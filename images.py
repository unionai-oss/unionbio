import re
from pathlib import Path
from flytekit import ImageSpec
from flytekit.image_spec.image_spec import ImageBuildEngine
from unionbio.config import current_registry

project_rt = Path(__file__).parent
prod_rt = project_rt.joinpath("src")
ws_rt = project_rt.joinpath("workspaces")

union_version = "union==0.1.85"

main_img = ImageSpec(
    name="main",
    platform="linux/amd64",
    python_version="3.12",
    packages=[union_version],
    source_root=prod_rt,
    conda_channels=["bioconda"],
    conda_packages=[
        "samtools",
        "bcftools",
        "bwa",
        "fastp",
        "hisat2",
        "bowtie2",
        "gatk4",
        "fastqc",
        "htslib",
        "multiqc",
    ],
    builder="fast-builder",
    registry=current_registry,
)

folding_img = ImageSpec(
    name="folding",
    platform="linux/amd64",
    python_version="3.12",
    packages=["transformers", "torch", union_version],
    source_root=prod_rt,
    conda_channels=["bioconda", "conda-forge"],
    conda_packages=[
        "prodigal",
        "biotite",
        "biopython",
        "py3Dmol",
        "matplotlib",
    ],
    builder="fast-builder",
    registry=current_registry,
)

parabricks_img = ImageSpec(
    name="parabricks",
    base_image="nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1",
    platform="linux/amd64",
    python_version="3.12",
    packages=[union_version],
    source_root=prod_rt,
    builder="fast-builder",
    registry=current_registry,
)

colabfold_img = ImageSpec(
    name="colabfold",
    platform="linux/amd64",
    apt_packages=["curl", "tar", "zstd", "gpg", "git"],
    python_version="3.10",
    packages=[
        union_version,
        # "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold",
        "colabfold[alphafold]",
        "jax[cuda12]",
        "tensorflow",
        "silence_tensorflow",
        "flytekitplugins-pod",
        "graphein",
    ],
    source_root=prod_rt,
    conda_channels=["bioconda", "conda-forge"],
    conda_packages=[
        "openmm",
        "pdbfixer",
        "kalign2",
        "hhsuite",
        "mmseqs2",
    ],
    commands=[
        # Install gcloud
        'echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] '
        'https://packages.cloud.google.com/apt cloud-sdk main" | tee -a '
        "/etc/apt/sources.list.d/google-cloud-sdk.list && curl "
        "https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor "
        "-o /usr/share/keyrings/cloud.google.gpg && apt-get update -y && "
        "apt-get install google-cloud-cli -y"
    ],
    builder="fast-builder",
    registry=current_registry,
)

build_scope = [
    # "main_img",
    # "folding_img",
    # "parabricks_img",
    "colabfold_img",
]


def update_img_config(config_path: Path, fqns: dict[str, str]):
    with open(config_path, "r") as f:
        cfg_content = f.read()

    for var, tag in fqns.items():
        cfg_content = re.sub(rf"{var} = .+", f'{var} = "{tag}"', cfg_content)

    with open(config_path, "w") as f:
        f.write(cfg_content)


def update_ws_config(config_path: Path, fqn: str):
    with open(config_path, "r") as file:
        yaml_data = file.readlines()

    # Manually parsing instead of using yaml library to preserve structure
    lines_out = []
    for line in yaml_data:
        if line.startswith("container_image:"):
            ll = line.split(": ")
            ll[1] = f"{fqn}\n"
            lines_out.append(": ".join(ll))
        else:
            lines_out.append(line)

    with open(config_path, "w") as file:
        file.writelines(lines_out)


## Entrypoints ##


def build():
    build_specs = []
    fqns = {}

    # Prepare builds
    for img_str in build_scope:
        spec = eval(img_str)
        fqns[f"{img_str}_fqn"] = spec.image_name()
        build_specs.append(spec)

    update_img_config(Path("src/unionbio/config.py"), fqns)

    for spec in build_specs:
        ImageBuildEngine().build(spec)


def build_test():
    build_specs = []
    fqns = {}

    # Prepare builds
    for img_str in build_scope:
        spec = eval(img_str)
        spec.tag_format = "{spec_hash}-test"
        spec.source_root = project_rt
        ns = spec.with_packages(["pytest"]).with_apt_packages(["screen", "htop"])
        fqn = ns.image_name()
        fqns[f"{img_str}_test_fqn"] = fqn
        build_specs.append(ns)

        # Update workspace config
        try:
            ws_cfg_path = ws_rt.joinpath(img_str.replace("_img", ""), "ws.yaml")
            update_ws_config(ws_cfg_path, fqn)
        except FileNotFoundError:
            print(f"Workspace config not found for {img_str}")

    update_img_config(Path("tests/config.py"), fqns)

    for spec in build_specs:
        ImageBuildEngine().build(spec)

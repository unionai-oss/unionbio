import re
import os
from pathlib import Path
from flytekit import ImageSpec
from flytekit.image_spec.image_spec import ImageBuildEngine

current_registry = os.getenv("IMAGE_SPEC_REGISTRY", "docker.io/unionbio")

project_rt = Path(__file__).parent
prod_rt = project_rt.joinpath("src")
ws_rt = project_rt.joinpath("workspaces")

union_version = "union==0.1.73"

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
    apt_packages=["curl"],
    builder="fast-builder",
    registry=current_registry,
)

alphafold_img = ImageSpec(
    name="alphafold",
    base_image="docker.io/unionbio/alphafold:base-20240910",
    platform="linux/amd64",
    python_version="3.12",
    packages=[union_version],
    apt_packages=["aria2"],
    source_root=prod_rt,
    entrypoint=[],
    # env={"PYTHONPATH": "/root:/opt/conda/lib/python3.11/site-packages"}, # Enable package discovery from base image
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

build_scope = [
    # "main_img",
    # "folding_img",
    "alphafold_img",
    # "parabricks_img",
]


def update_img_config(config_path: Path, fqns: dict[str, str]):
    with open(config_path, "r") as f:
        cfg_content = f.read()

    for var, tag in fqns.items():
        cfg_content = re.sub(rf"{var} = .+", f'{var} = "{tag}"', cfg_content)

    with open(config_path, "w") as f:
        f.write(cfg_content)


def update_ws_config(config_path: Path, fqn: str):
    with open(config_path, 'r') as file:
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

    with open(config_path, 'w') as file:
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
        spec.packages.append("pytest")
        fqn = spec.image_name()
        fqns[f"{img_str}_test_fqn"] = fqn
        build_specs.append(spec)

        # Update workspace config
        try:
            ws_cfg_path = ws_rt.joinpath(img_str.replace("_img", ""), "ws.yaml")
            update_ws_config(ws_cfg_path, fqn)
        except FileNotFoundError:
            print(f"Workspace config not found for {img_str}")

    update_img_config(Path("tests/config.py"), fqns)

    for spec in build_specs:
        ImageBuildEngine().build(spec)

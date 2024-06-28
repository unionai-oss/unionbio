import re
from pathlib import Path
from flytekit import ImageSpec
from flytekit.image_spec.image_spec import ImageBuildEngine

current_registry = "docker.io/unionbio"
test_rt = Path(__file__).parent
prod_rt = test_rt.joinpath("src")

main_img = ImageSpec(
    name="main",
    platform="linux/amd64",
    python_version="3.11",
    packages=["flytekit", "unionai"],
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
    python_version="3.11",
    packages=["flytekit", "transformers", "torch", "unionai"],
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
    python_version="3.10",
    packages=["flytekit", "unionai"],
    source_root=prod_rt,
    builder="fast-builder",
    registry=current_registry,
)

build_scope = [
    "main_img",
    # "folding_img",
    # "parabricks_img",
]


def update_img_config(config_path: Path, fqns: dict[str, str]):
    with open(config_path, "r") as f:
        cfg_content = f.read()

    for var, tag in fqns.items():
        cfg_content = re.sub(rf"{var} = .+", f'{var} = "{tag}"', cfg_content)

    with open(config_path, "w") as f:
        f.write(cfg_content)

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
        spec.name = f"{spec.name}-test"
        spec.source_root = test_rt
        spec.packages.append("pytest")
        fqns[f"{img_str}_test_fqn"] = spec.image_name()
        build_specs.append(spec)

    update_img_config(Path("tests/config.py"), fqns)

    for spec in build_specs:
        ImageBuildEngine().build(spec)

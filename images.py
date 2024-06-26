import re
from pathlib import Path
from flytekit import ImageSpec
from flytekit.image_spec.image_spec import ImageBuildEngine

current_registry = "localhost:30000"
test_rt = Path(__file__).parent
prod_rt = test_rt.joinpath("src")

main_img = ImageSpec(
    name="unionbio-main",
    source_root=prod_rt,
    packages=["flytekit"],
    python_version="3.11",
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
        "pytest", # Temporary
    ],
    builder="fast-builder",
    registry=current_registry,
)

folding_img = ImageSpec(
    name="unionbio-protein",
    platform="linux/amd64",
    python_version="3.11",
    source_root=prod_rt,
    packages=["flytekit", "transformers", "torch"],
    conda_channels=["bioconda", "conda-forge"],
    conda_packages=[
        "prodigal",
        "biotite",
        "biopython",
        "py3Dmol",
        "matplotlib",
    ],
    registry=current_registry,
)

parabricks_img = ImageSpec(
    name="unionbio-parabricks",
    base_image="nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1",
    source_root=prod_rt,
    python_version="3.10",
    packages=["flytekit"],
    registry=current_registry,
)

build_scope = [
    "main_img",
    "folding_img",
    "parabricks_img",
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

import os
from flytekit import ImageSpec
from flytekit.image_spec.image_spec import ImageBuildEngine
from unionbio.config import prod_rt

union_version = "union==0.1.103"
current_registry = os.getenv("IMAGE_SPEC_REGISTRY", "docker.io/unionbio")

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
    registry=current_registry,
)

parabricks_img = ImageSpec(
    name="parabricks",
    base_image="nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1",
    platform="linux/amd64",
    python_version="3.12",
    packages=[union_version],
    source_root=prod_rt,
    registry=current_registry,
)

colabfold_img = ImageSpec(
    name="colabfold",
    platform="linux/amd64",
    apt_packages=[
        "curl",
        "tar",
        "zstd",
        "gpg",
        "git",
        "gcc",
        "wget",
        "unzip",
        "build-essential",
        "libc6-dev",
    ],
    python_version="3.10",
    packages=[
        union_version,
        "flytekitplugins-pod",
        "graphein",
        "zstandard",
        "colabfold[alphafold-minus-jax]@git+https://github.com/sokrypton/ColabFold.git",
        "colabfold[alphafold]",
        "jax[cuda12]==0.4.35",
    ],
    conda_channels=["conda-forge", "bioconda"],
    conda_packages=[
        "openmm==7.7.0",
        "pdbfixer",
        "kalign2==2.04",
        "hhsuite==3.3.0",
        "mmseqs2==15.6f452",
    ],
    source_root=prod_rt,
    commands=[
        # Install gcloud
        'echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
        curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
        apt-get update -y && apt-get install google-cloud-cli -y',
    ],
    registry=current_registry,
)

# Determines which images will be built for prod or test
build_scope = [
    "main_img",
    "parabricks_img",
    "colabfold_img",
]


def build():
    # Prepare builds
    for img_str in build_scope:
        spec = eval(img_str)
        ImageBuildEngine().build(spec)

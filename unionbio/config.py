import logging
from typing import List, Dict
from pprint import pprint
from pathlib import Path
from flytekit import ImageSpec

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(
    logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s")
)
# logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

folding_img = ImageSpec(
    name="unionbio-protein",
    base_image="ghcr.io/flyteorg/flytekit:py3.11-1.12.0",
    packages=["requests"],
    # conda_channels=["bioconda"],
    # conda_packages=["biopython", "biotite", "py3Dmol"],
    registry="ghcr.io/unionai-oss"
)

parabricks_img = ImageSpec(
    name="unionbio-parabricks",
    base_image="nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1",
    python_version="3.10",
    packages=["flytekit"],
    registry="ghcr.io/unionai-oss"
)

# main_img = ImageSpec(
#     name="unionbio-main",
#     base_image="ghcr.io/flyteorg/flytekit:py3.11-1.12.0",
#     # conda_channels=["bioconda"],
#     conda_packages=["requests"],
#     registry="ghcr.io/unionai-oss",
#     platform="linux/amd64"
# )

main_img = ImageSpec(
    base_image="ubuntu:20.04",
    python_version="3.11",
    packages=["flytekit"],
    conda_packages=["pytorch", "cpuonly"],
    conda_channels=["pytorch"],
)

# Define parabricks image
pb_image = "ghcr.io/unionai/dgx-parabricks:20240416"

seq_dir_pth = "s3://my-s3-bucket/my-data/sequences"
ref_loc = "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
ref_hash = str(hash(ref_loc))[:4]

# Tool config
fastp_cpu = "3"


import logging
import tomllib
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

env_path = Path(__file__).parent.parent.joinpath("pixi.toml")
env = {}
with open(env_path, "rb") as f:
    env = tomllib.load(f)

def make_paks(deps: Dict) -> List[str]:
    conda_paks = []
    for k,v in deps.items():
        if v == "*":
            conda_paks.append(k)
        else:
            conda_paks.append(f"{k}{v}")
    return conda_paks

# Define main image
main_image = ImageSpec(
    name="unionbio-main",
    base_image="ghcr.io/flyteorg/flytekit:latest",
    conda_channels=env['feature']['main']['channels'],
    conda_packages=make_paks(env['feature']['main']['dependencies']),
    registry="ghcr.io/unionai-oss",
)

# Define parabricks image
pb_image = "ghcr.io/unionai/dgx-parabricks:20240416"

seq_dir_pth = "s3://my-s3-bucket/my-data/sequences"
ref_loc = "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
ref_hash = str(hash(ref_loc))[:4]

# Tool config
fastp_cpu = "3"


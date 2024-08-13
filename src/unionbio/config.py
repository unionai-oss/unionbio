import os
import logging
from pathlib import Path

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(
    logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s")
)
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

# Default paths
remote_ref = ''
remote_reads = ''
remote_sites = ''
ref_hash = str(hash(remote_ref))[:4]

# Tool config
fastp_cpu = "3"

current_registry = os.getenv("IMAGE_SPEC_REGISTRY", "docker.io/unionbio")
src_rt = Path(__file__).parent.parent

# Image tags
# While tasks can reference imageSpec directly, using the tag allows registering tasks
# from a containerized environment. These also contain the actual unionbio package.
main_img_fqn = "docker.io/unionbio/main:NqPVgtIraS7qOGkEFbTjWQ"
folding_img_fqn = "docker.io/unionbio/folding:kDb3vJxS4f_uisEsy833SA"
parabricks_img_fqn = "docker.io/unionbio/parabricks:d3RaRIVpChR_qkAbDBIylg"

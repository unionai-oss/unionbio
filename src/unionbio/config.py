import logging
from pathlib import Path

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(
    logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s")
)
# logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

seq_dir_pth = "s3://my-s3-bucket/my-data/sequences"
ref_loc = "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
ref_hash = str(hash(ref_loc))[:4]

# Tool config
fastp_cpu = "3"

# current_registry = "ghcr.io/unionai-oss"
current_registry = "localhost:30000"
src_rt = Path(__file__).parent.parent

# Image tags
# While tasks can reference imageSpec directly, using the tag allows registering tasks
# from a containerized environment. These also contain the actual unionbio package.
main_img_fqn = "localhost:30000/unionbio-main:rpC35uIsi7b7e_EZh9atRA"
folding_img_fqn = "localhost:30000/unionbio-protein:LBzgRjOPHBOLPB6In5X_sA"
parabricks_img_fqn = "localhost:30000/unionbio-parabricks:NRCm_E3moyA0Li29GlR8Nw"

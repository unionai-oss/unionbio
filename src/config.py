import logging

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(
    logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s")
)
# logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

base_image = "localhost:30000/variant-discovery:latest"
seq_dir_pth = "s3://my-s3-bucket/my-data/sequences"
ref_loc = "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
ref_hash = str(hash(ref_loc))[:4]

test_assets = {
    'local_seq_dir': '/mnt/sequences',
    'local_tiny_seq_dir': '/mnt/tiny_sequences',
    'local_filt_dir': '/mnt/filtered_sequences',
    'local_sam_dir': '/mnt/alignments',
}

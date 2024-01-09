import logging

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(
    logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s")
)
# logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

base_image = "ghcr.io/pryce-turner/variant-discovery:20240102"
seq_dir_pth = "s3://my-s3-bucket/my-data/sequences"
ref_loc = "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
ref_hash = str(hash(ref_loc))[:4]

test_assets = {
    "seq_dir": "/root/src/tests/assets/sequences",
    "filt_dir": "/root/src/tests/assets/filtered_sequences",
    "bt2_sam_dir": "/root/src/tests/assets/alignments/bt2",
    "hs2_sam_dir": "/root/src/tests/assets/alignments/hs2",
    "ref_path": "/root/src/tests/assets/references/GRCh38_short.fasta",
    "idx_dir": "/root/src/tests/assets/indices",
    "bt2_idx_dir": "/root/src/tests/assets/indices/bt2",
    "hs2_idx_dir": "/root/src/tests/assets/indices/hs2",
    "fastqc_dir": "/root/src/tests/assets/fastqc",
}

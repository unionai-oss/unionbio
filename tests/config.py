from pathlib import Path

proj_rt = Path(__file__).parent.parent
test_dir = proj_rt.joinpath("tests")

test_assets = {
    "raw_seq_dir": f"{test_dir}/assets/sequences/raw",
    "raw_read": f"{test_dir}/assets/sequences/raw/ERR250683-tiny_1.fastq.gz",
    "filt_seq_dir": f"{test_dir}/assets/sequences/filtered",
    "folding_seq_dir": f"{test_dir}/assets/sequences/folding",
    "bt2_sam_dir": f"{test_dir}/assets/alignments/bt2",
    "hs2_sam_dir": f"{test_dir}/assets/alignments/hs2",
    "pb_fq2bam_dir": f"{test_dir}/assets/alignments/pb_fq2bam",
    "pb_haplocall_dir": f"{test_dir}/assets/alignments/pb_haplocall",
    "sort_dir": f"{test_dir}/assets/alignments/sorted",
    "dedup_dir": f"{test_dir}/assets/alignments/deduped",
    "ref_path": f"{test_dir}/assets/references/GRCh38_short.fasta",
    "ref_idx_path": f"{test_dir}/assets/references/GRCh38_short.fasta.fai",
    "ref_dir": f"{test_dir}/assets/references/",
    "ref_fn": "GRCh38_short.fasta",
    "sites_path": f"{test_dir}/assets/sites/known_indels_trunc.hg38.vcf.gz",
    "sites_idx_path": f"{test_dir}/assets/sites/known_indels_trunc.hg38.vcf.gz.tbi",
    "vcf_path": f"{test_dir}/assets/vcfs/test-sample-1_test-caller.vcf.gz",
    "vcf_idx_path": f"{test_dir}/assets/vcfs/test-sample-1_test-caller.vcf.gz.tbi",
    "vcf_dir": f"{test_dir}/assets/vcfs/",
    "idx_dir": f"{test_dir}/assets/indices",
    "bt2_idx_dir": f"{test_dir}/assets/indices/bt2",
    "hs2_idx_dir": f"{test_dir}/assets/indices/hs2",
    "bwa_idx_dir": f"{test_dir}/assets/indices/bwa",
    "fastqc_dir": f"{test_dir}/assets/fastqc",
    "prot_path": f"{test_dir}/assets/proteins/folding_proteins.fasta",
    "genes_path": f"{test_dir}/assets/genes/folding_genes.gff",
}

# Image tags
# While tasks can reference imageSpec directly, using the tag allows registering tasks
# from a containerized environment. These also contain the actual unionbio package.
main_img_test_fqn = "docker.io/unionbio/main:vpYJ9uAzyUlqo3BO4_VoLA-test"
folding_img_test_fqn = "docker.io/unionbio/folding:cHtEowX2F_skDPfJb4dvaA-test"
parabricks_img_test_fqn = "docker.io/unionbio/parabricks:Y6cllYwwhvHQgZOn7zClQw-test"

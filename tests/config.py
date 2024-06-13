from pathlib import Path

proj_rt = Path(__file__).parent.parent
test_dir = proj_rt.joinpath("tests")

test_assets = {
    "raw_seq_dir": f"{test_dir}/assets/sequences/raw",
    "filt_seq_dir": f"{test_dir}/assets/sequences/filtered",
    "folding_seq_dir": f"{test_dir}/assets/sequences/folding",
    "bt2_sam_dir": f"{test_dir}/assets/alignments/bt2",
    "hs2_sam_dir": f"{test_dir}/assets/alignments/hs2",
    "pb_fq2bam_dir": f"{test_dir}/assets/alignments/pb_fq2bam",
    "pb_haplocall_dir": f"{test_dir}/assets/alignments/pb_haplocall",
    "sort_dir": f"{test_dir}/assets/alignments/sorted",
    "dedup_dir": f"{test_dir}/assets/alignments/deduped",
    "ref_path": f"{test_dir}/assets/references/GRCh38_short.fasta",
    "ref_dir": f"{test_dir}/assets/references/",
    "ref_fn": "GRCh38_short.fasta",
    "sites_path": f"{test_dir}/assets/sites/known_indels_trunc.hg38.vcf.gz",
    "sites_idx_path": f"{test_dir}/assets/sites/known_indels_trunc.hg38.vcf.gz.tbi",
    "vcf_path": f"{test_dir}/assets/vcfs/test-sample_test-caller.vcf.gz",
    "vcf_idx_path": f"{test_dir}/assets/vcfs/test-sample_test-caller.vcf.gz.tbi",
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
main_img_test_fqn = "localhost:30000/unionbio-main-test:z_VQI0cxFSK_kj7GgYpbsQ"
folding_img_test_fqn = "localhost:30000/unionbio-protein-test:NlNhyap5cebT0lMp_D4hmQ"
parabricks_img_test_fqn = "localhost:30000/unionbio-parabricks-test:tjWURFYPRxGXDlo8r09GKQ"
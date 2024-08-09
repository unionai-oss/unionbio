from pathlib import Path

proj_rt = Path(__file__).parent.parent
test_dir = proj_rt.joinpath("tests")

test_assets = {
    "raw_seq_dir": f"{test_dir}/assets/sequences/raw",
    "raw_read": f"{test_dir}/assets/sequences/raw/ERR250683-tiny_1.fastq.gz",
    "filt_seq_dir": f"{test_dir}/assets/sequences/filtered",
    "folding_seq_dir": f"{test_dir}/assets/sequences/folding",
    "al_dir": f"{test_dir}/assets/alignments",
    "bt2_sam_dir": f"{test_dir}/assets/alignments/bt2",
    "hs2_sam_dir": f"{test_dir}/assets/alignments/hs2",
    "bwa_sam_path": f"{test_dir}/assets/alignments/bwa/ERR250683-tiny_bwa.sam",
    "pb_fq2bam_dir": f"{test_dir}/assets/alignments/pb_fq2bam",
    "pb_haplocall_dir": f"{test_dir}/assets/alignments/pb_haplocall",
    "sort_dir": f"{test_dir}/assets/alignments/sorted",
    "bam_dir": f"{test_dir}/assets/alignments/recal/bam",
    "dedup_dir": f"{test_dir}/assets/alignments/deduped",
    "recal_in": f"{test_dir}/assets/alignments/deduped/chr21",
    "recal_out": f"{test_dir}/assets/alignments/recal/SRR812824-chr21_bowtie2_sorted_deduped_recal.sam",
    "ref_path": f"{test_dir}/assets/references/GRCh38_short.fasta",
    "ref_idx_path": f"{test_dir}/assets/references/GRCh38_short.fasta.fai",
    "ref_dir": f"{test_dir}/assets/references/",
    "ref_fn": "GRCh38_short.fasta",
    "sites_dir": f"{test_dir}/assets/sites",
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
main_img_test_fqn = "docker.io/unionbio/main:4wKChuc4z1UC5Y5CljlzTg-test"
folding_img_test_fqn = "docker.io/unionbio/folding:ZzkXAQAcVudjddDAymIBgg-test"
parabricks_img_test_fqn = "docker.io/unionbio/parabricks:PvCm4VJqXalHgP5_aII0kw-test"

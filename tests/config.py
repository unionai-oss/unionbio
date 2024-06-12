from pathlib import Path

cwd = Path(__file__).parent

test_assets = {
    "raw_seq_dir": f"{cwd}/assets/sequences/raw",
    "filt_seq_dir": f"{cwd}/assets/sequences/filtered",
    "folding_seq_dir": f"{cwd}/assets/sequences/folding",
    "bt2_sam_dir": f"{cwd}/assets/alignments/bt2",
    "hs2_sam_dir": f"{cwd}/assets/alignments/hs2",
    "pb_fq2bam_dir": f"{cwd}/assets/alignments/pb_fq2bam",
    "pb_haplocall_dir": f"{cwd}/assets/alignments/pb_haplocall",
    "sort_dir": f"{cwd}/assets/alignments/sorted",
    "dedup_dir": f"{cwd}/assets/alignments/deduped",
    "ref_path": f"{cwd}/assets/references/GRCh38_short.fasta",
    "ref_dir": f"{cwd}/assets/references/",
    "ref_fn": "GRCh38_short.fasta",
    "sites_path": f"{cwd}/assets/sites/known_indels_trunc.hg38.vcf.gz",
    "sites_idx_path": f"{cwd}/assets/sites/known_indels_trunc.hg38.vcf.gz.tbi",
    "vcf_path": f"{cwd}/assets/vcfs/test-sample_test-caller.vcf.gz",
    "vcf_idx_path": f"{cwd}/assets/vcfs/test-sample_test-caller.vcf.gz.tbi",
    "vcf_dir": f"{cwd}/assets/vcfs/",
    "idx_dir": f"{cwd}/assets/indices",
    "bt2_idx_dir": f"{cwd}/assets/indices/bt2",
    "hs2_idx_dir": f"{cwd}/assets/indices/hs2",
    "bwa_idx_dir": f"{cwd}/assets/indices/bwa",
    "fastqc_dir": f"{cwd}/assets/fastqc",
    "prot_path": f"{cwd}/assets/proteins/folding_proteins.fasta",
    "genes_path": f"{cwd}/assets/genes/folding_genes.gff",
}

# Image tags
# While tasks can reference imageSpec directly, using the tag allows registering tasks
# from a containerized environment. These also contain the actual unionbio package.
main_img_test_fqn = "localhost:30000/unionbio-main-test:JEgls19QqlpRRUSV1AoG_w"
folding_img_test_fqn = "localhost:30000/unionbio-protein-test:4OwU32HTc9VEYteE4vxhbw"
parabricks_img_test_fqn = "localhost:30000/unionbio-parabricks-test:61aLsbFRwZsvd0JvIXRTEA"
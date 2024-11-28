from pathlib import Path

proj_rt = Path(__file__).parent.parent
test_dir = proj_rt.joinpath("tests")
assets = test_dir.joinpath("assets")

test_assets = {
    "raw_seq_dir": assets.joinpath("sequences/raw"),
    "filt_seq_dir": assets.joinpath("sequences/filtered"),
    "folding_seq_dir": assets.joinpath("sequences/folding"),
    "al_dir": assets.joinpath("alignments"),
    "bt2_sam_dir": assets.joinpath("alignments/bt2"),
    "hs2_sam_dir": assets.joinpath("alignments/hs2"),
    "pb_fq2bam_dir": assets.joinpath("alignments/pb_fq2bam"),
    "pb_haplocall_dir": assets.joinpath("alignments/pb_haplocall"),
    "sort_dir": assets.joinpath("alignments/sorted"),
    "dedup_dir": assets.joinpath("alignments/deduped"),
    "bam_dir": assets.joinpath("alignments/bam"),
    "ref_dir": assets.joinpath("references"),
    "sites_dir": assets.joinpath("sites"),
    "vcf_dir": assets.joinpath("vcfs"),
    "isec_dir": assets.joinpath("vcfs/isec"),
    "idx_dir": assets.joinpath("indices"),
    "bt2_idx_dir": assets.joinpath("indices/bt2"),
    "hs2_idx_dir": assets.joinpath("indices/hs2"),
    "bwa_idx_dir": assets.joinpath("indices/bwa"),
    "fastqc_dir": assets.joinpath("fastqc"),
    "ref_fn": "GRCh38_short.fasta",
    "prot_path": assets.joinpath("proteins/folding_proteins.fasta"),
    "ref_path": assets.joinpath("references/GRCh38_short.fasta"),
    "ref_idx_path": assets.joinpath("references/GRCh38_short.fasta.fai"),
    "bwa_sam_path": assets.joinpath("alignments/bwa/ERR250683-tiny_bwa.sam"),
    "recal_out": assets.joinpath(
        "alignments/recal/SRR812824-chr21_bowtie2_sorted_deduped_recal.sam"
    ),
    "sites_path": assets.joinpath("sites/known_indels_trunc.hg38.vcf.gz"),
    "sites_idx_path": assets.joinpath("sites/known_indels_trunc.hg38.vcf.gz.tbi"),
    "vcf_path": assets.joinpath("vcfs/test-sample-1_test-caller.vcf.gz"),
    "vcf_idx_path": assets.joinpath("vcfs/test-sample-1_test-caller.vcf.gz.tbi"),
    "raw_read": assets.joinpath("sequences/raw/ERR250683-tiny_1.fastq.gz"),
    "protein_path": assets.joinpath("proteins"),
}

# Image tags
# While tasks can reference imageSpec directly, using the tag allows registering tasks
# from a containerized environment. These also contain the actual unionbio package.
main_img_test_fqn = "docker.io/unionbio/main:jnnHUFNpX1D6ogVmahU2FQ-test"
parabricks_img_test_fqn = "docker.io/unionbio/parabricks:HR03xNoQHQVS3RAGwiMBWw-test"
colabfold_img_test_fqn = "docker.io/unionbio/colabfold:_Q52hZrWIN9CbJgX4Yo6ow-test"

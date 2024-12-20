import os
from pathlib import Path
from tests.config import test_assets
from tests.utils import copy_dir_conts
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from unionbio.types import Alignment, Reads, Reference, VCF, Protein


def test_raw_sample_fname():
    o1, o2 = Reads("test").get_read_fnames()
    assert o1 == "test_1.fastq.gz"


def test_raw_sample_make_all():
    samps = Reads.make_all(test_assets["raw_seq_dir"])
    assert len(samps) == 1
    assert isinstance(samps[0], Reads)
    assert "ERR250683-tiny_1.fastq.gz" in samps[0].read1.path
    assert "ERR250683-tiny_2.fastq.gz" in samps[0].read2.path


def test_filt_sample_fname():
    filt_samp = Reads.make_all(test_assets["filt_seq_dir"])[0]
    o1, o2 = filt_samp.get_read_fnames()
    rep = filt_samp.get_report_fname()
    assert o1 == "ERR250683-tiny_1.filt.fastq.gz"
    assert rep == "ERR250683-tiny_fastq-filter-report.json"


def test_filt_sample_make_all():
    filt_samps = Reads.make_all(test_assets["filt_seq_dir"])
    assert len(filt_samps) == 1
    assert isinstance(filt_samps[0], Reads)
    assert filt_samps[0].sample == "ERR250683-tiny"
    assert "ERR250683-tiny_1.filt.fastq.gz" in filt_samps[0].read1.path
    assert "ERR250683-tiny_2.filt.fastq.gz" in filt_samps[0].read2.path
    assert "ERR250683-tiny_fastq-filter-report.json" in filt_samps[0].filt_report.path


def test_reads_aggregate(tmp_path):
    copy_dir_conts(test_assets["raw_seq_dir"], tmp_path)
    samp = Reads.make_all(Path(tmp_path))[0]
    target = samp.aggregate(target=Path("/tmp/some/other/dir"))
    assert Path(target).exists()
    assert samp.read2.path == Path(target).joinpath("ERR250683-tiny_2.fastq.gz")
    assert samp.read1.path == Path(target).joinpath("ERR250683-tiny_1.fastq.gz")


def test_alignment_file_fname():
    al = Alignment("test", "bowtie2", "sam").get_alignment_fname()
    assert al == "test_bowtie2_aligned.sam"


def test_alignment_file_make_all():
    sams = Alignment.make_all(test_assets["bt2_sam_dir"])
    assert len(sams) == 1
    assert isinstance(sams[0], Alignment)
    assert sams[0].sample == "ERR250683-tiny"
    assert sams[0].aligner in ["hisat2", "bowtie2"]
    assert "ERR250683-tiny_bowtie2_aligned.sam" in sams[0].alignment.path
    assert "ERR250683-tiny_bowtie2_aligned_report.txt" in sams[0].alignment_report.path


def test_alignment_aggregate(tmp_path):
    copy_dir_conts(test_assets["bt2_sam_dir"], tmp_path)
    al = Alignment.make_all(Path(tmp_path))[0]
    target = al.aggregate(target=Path("/tmp/some/other/dir"))
    assert Path(target).exists()
    assert al.alignment.path == Path(target).joinpath(
        "ERR250683-tiny_bowtie2_aligned.sam"
    )
    assert al.alignment_report.path == Path(target).joinpath(
        "ERR250683-tiny_bowtie2_aligned_report.txt"
    )


def test_reference():
    ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=test_assets["ref_dir"]))
    assert isinstance(ref.ref_dir, FlyteDirectory)
    assert ref.ref_name in os.listdir(ref.ref_dir.path)
    assert ref.get_ref_path() == test_assets["ref_dir"].joinpath(test_assets["ref_fn"])


def test_remote_reference():
    ref = Reference.from_remote(
        url="https://raw.githubusercontent.com/unionai-oss/unionbio/main/tests/assets/references/GRCh38_short.fasta"
    )
    assert isinstance(ref, Reference)
    assert ref.ref_name == "GRCh38_short.fasta"
    assert Path(ref.ref_dir.path).exists()


def test_reference_aggregate(tmp_path):
    tp1 = Path(tmp_path).joinpath("d1")
    tp2 = Path(tmp_path).joinpath("d2")
    copy_dir_conts(test_assets["ref_dir"], tp1)
    ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=tp1))
    target = ref.aggregate(target=tp2)
    assert Path(target).exists()
    assert ref.ref_dir.path == target
    assert ref.ref_name in os.listdir(target)


def test_vcf():
    vcf = VCF(
        "test-sample",
        "test-caller",
        FlyteFile(path=test_assets["vcf_path"]),
        FlyteFile(path=test_assets["vcf_idx_path"]),
    )
    assert vcf.sample == "test-sample"
    assert vcf.caller == "test-caller"
    assert vcf.vcf.path == test_assets["vcf_path"]
    assert vcf.vcf_idx.path == test_assets["vcf_idx_path"]
    assert vcf.get_vcf_fname() == "test-sample_test-caller.g.vcf.gz"
    assert vcf.get_vcf_idx_fname() == "test-sample_test-caller.g.vcf.gz.tbi"


def test_vcf_make_all():
    vcfs = VCF.make_all(test_assets["vcf_dir"], include=["*SRR812824-chr21*"])
    assert len(vcfs) == 1
    assert isinstance(vcfs[0], VCF)
    assert vcfs[0].sample == "SRR812824-chr21"
    assert vcfs[0].caller == "haplo-caller"
    assert vcfs[0].vcf.path == str(test_assets["vcf_path"])
    assert vcfs[0].vcf_idx.path == str(test_assets["vcf_idx_path"])
    assert vcfs[0].get_vcf_fname() == "SRR812824-chr21_haplo-caller.g.vcf.gz"
    assert vcfs[0].get_vcf_idx_fname() == "SRR812824-chr21_haplo-caller.g.vcf.gz.tbi"


def test_protein():
    prot = Protein("my-monomer", FlyteFile(path="sequence-path"))
    assert prot.sample == "my-monomer"
    assert prot.sequence.path == "sequence-path"
    assert prot.get_prot_fname() == "my-monomer.fasta"
    assert prot.get_genes_fname() == "my-monomer.gff"


def test_protein_make_all(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    p = Protein.make_all(tmp_path)[0]
    assert isinstance(p, Protein)

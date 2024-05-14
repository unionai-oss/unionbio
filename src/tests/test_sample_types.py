import os
from pathlib import Path
from config import test_assets
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from datatypes.alignment import Alignment
from datatypes.reads import Reads
from datatypes.reference import Reference
from datatypes.variants import VCF
from datatypes.protein import Protein


def test_raw_sample_fname():
    o1, o2 = Reads("test").get_read_fnames()
    assert o1 == "test_1.fastq.gz"


def test_raw_sample_make_all():
    samps = Reads.make_all(Path(test_assets["raw_seq_dir"]))
    assert len(samps) == 1
    assert isinstance(samps[0], Reads)
    assert "ERR250683-tiny_1.fastq.gz" in samps[0].read1.path
    assert "ERR250683-tiny_2.fastq.gz" in samps[0].read2.path


def test_filt_sample_fname():
    filt_samp = Reads.make_all(Path(test_assets["filt_seq_dir"]))[0]
    o1, o2 = filt_samp.get_read_fnames()
    rep = filt_samp.get_report_fname()
    assert o1 == "ERR250683-tiny_1.filt.fastq.gz"
    assert rep == "ERR250683-tiny_fastq-filter-report.json"


def test_filt_sample_make_all():
    filt_samps = Reads.make_all(Path(test_assets["filt_seq_dir"]))
    assert len(filt_samps) == 1
    assert isinstance(filt_samps[0], Reads)
    assert filt_samps[0].sample == "ERR250683-tiny"
    assert "ERR250683-tiny_1.filt.fastq.gz" in filt_samps[0].read1.path
    assert "ERR250683-tiny_2.filt.fastq.gz" in filt_samps[0].read2.path
    assert "ERR250683-tiny_fastq-filter-report.json" in filt_samps[0].filt_report.path


def test_alignment_file_fname():
    al = Alignment("test", "bowtie2", "sam").get_alignment_fname()
    assert al == "test_bowtie2_aligned.sam"


def test_alignment_file_make_all():
    sams = Alignment.make_all(Path(test_assets["bt2_sam_dir"]))
    assert len(sams) == 1
    assert isinstance(sams[0], Alignment)
    assert sams[0].sample == "ERR250683-tiny"
    assert sams[0].aligner in ["hisat2", "bowtie2"]
    assert "ERR250683-tiny_bowtie2_aligned.sam" in sams[0].alignment.path
    assert "ERR250683-tiny_bowtie2_aligned_report.txt" in sams[0].alignment_report.path


def test_reference():
    ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=test_assets["ref_dir"]))
    assert isinstance(ref.ref_dir, FlyteDirectory)
    assert ref.ref_name in os.listdir(ref.ref_dir.path)
    assert ref.get_ref_path() == Path(test_assets["ref_dir"]).joinpath(
        test_assets["ref_fn"]
    )


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
    assert vcf.get_vcf_fname() == "test-sample_test-caller.vcf.gz"
    assert vcf.get_vcf_idx_fname() == "test-sample_test-caller.vcf.gz.tbi"



def test_vcf_make_all():
    vcfs = VCF.make_all(Path(test_assets["vcf_dir"]))
    assert len(vcfs) == 1
    assert isinstance(vcfs[0], VCF)
    assert vcfs[0].sample == "test-sample"
    assert vcfs[0].caller == "test-caller"
    assert vcfs[0].vcf.path == test_assets["vcf_path"]
    assert vcfs[0].vcf_idx.path == test_assets["vcf_idx_path"]
    assert vcfs[0].get_vcf_fname() == "test-sample_test-caller.vcf.gz"
    assert vcfs[0].get_vcf_idx_fname() == "test-sample_test-caller.vcf.gz.tbi"


def test_protein():
    prot = Protein("test-protein", FlyteFile(path="test-path"))
    assert prot.name == "test-protein"
    assert prot.protein.path == "test-path"
    assert prot.get_prot_fname() == "test-protein_proteins.fasta"
    assert prot.get_genes_fname() == "test-protein_genes.gff"
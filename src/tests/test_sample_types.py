from pathlib import Path
from config import test_assets
from datatypes.alignment import Alignment
from datatypes.reads import Reads


def test_raw_sample_fname():
    o1, o2 = Reads("test").get_read_fnames()
    assert o1 == "test_1.fastq.gz"


def test_raw_sample_make_all():
    samps = Reads.make_all(Path(test_assets["seq_dir"]))
    assert len(samps) == 1
    assert isinstance(samps[0], Reads)
    assert "ERR250683-tiny_1.fastq.gz" in samps[0].read1.path
    assert "ERR250683-tiny_2.fastq.gz" in samps[0].read2.path


def test_filt_sample_fname():
    filt_samp = Reads.make_all(Path(test_assets["filt_dir"]))[0]
    o1, o2 = filt_samp.get_read_fnames()
    rep = filt_samp.get_report_fname()
    assert o1 == "ERR250683-tiny_1.filt.fastq.gz"
    assert rep == "ERR250683-tiny_fastq-filter-report.json"


def test_filt_sample_make_all():
    filt_samps = Reads.make_all(Path(test_assets["filt_dir"]))
    assert len(filt_samps) == 1
    assert isinstance(filt_samps[0], Reads)
    assert filt_samps[0].sample == "ERR250683-tiny"
    assert "ERR250683-tiny_1.filt.fastq.gz" in filt_samps[0].read1.path
    assert "ERR250683-tiny_2.filt.fastq.gz" in filt_samps[0].read2.path
    assert "ERR250683-tiny_fastq-filter-report.json" in filt_samps[0].filt_report.path


def test_sam_file_fname():
    sam = Alignment("test", "bowtie2").get_alignment_fname()
    assert sam == "test_bowtie2_aligned.sam"


def test_sam_file_make_all():
    sams = Alignment.make_all(Path(test_assets["bt2_sam_dir"]))
    assert len(sams) == 1
    assert isinstance(sams[0], Alignment)
    assert sams[0].sample == "ERR250683-tiny"
    assert sams[0].aligner in ["hisat2", "bowtie2"]
    assert "ERR250683-tiny_bowtie2_aligned.sam" in sams[0].sam.path
    assert "ERR250683-tiny_bowtie2_aligned_report.txt" in sams[0].report.path

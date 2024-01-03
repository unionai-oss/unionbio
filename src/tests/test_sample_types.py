from pathlib import Path
from config import test_assets
from tasks.sample_types import RawSample, FiltSample, SamFile


def test_raw_sample_fname():
    o1, o2 = RawSample("test").make_filenames()
    assert o1 == "test_1.fastq.gz"


def test_raw_sample_make_all():
    samps = RawSample.make_all(Path(test_assets["seq_dir"]))
    print(samps)
    assert len(samps) == 1
    assert isinstance(samps[0], RawSample)
    assert "ERR250683-tiny_1.fastq.gz" in samps[0].raw_r1.path
    assert "ERR250683-tiny_2.fastq.gz" in samps[0].raw_r2.path


def test_filt_sample_fname():
    o1, o2, rep = FiltSample("test").make_filenames()
    assert o1 == "test_1.filt.fastq.gz"


def test_filt_sample_make_all():
    filt_samps = FiltSample.make_all(Path(test_assets["filt_dir"]))
    assert len(filt_samps) == 1
    assert isinstance(filt_samps[0], FiltSample)
    assert filt_samps[0].sample == "ERR250683-tiny"
    assert "ERR250683-tiny_1.filt.fastq.gz" in filt_samps[0].filt_r1.path
    assert "ERR250683-tiny_2.filt.fastq.gz" in filt_samps[0].filt_r2.path
    assert "ERR250683-tiny_filt-report.json" in filt_samps[0].report.path


def test_sam_file_fname():
    o1, o2, rep = SamFile("test").make_filenames()
    assert o1 == "test_bowtie2_aligned.sam"


def test_sam_file_make_all():
    sams = SamFile.make_all(Path(test_assets["sam_dir"]))
    assert len(sams) == 1
    assert isinstance(sams[0], SamFile)
    assert sams[0].sample == "ERR250683-tiny"
    assert sams[0].aligner in ["hisat2", "bowtie2"]
    assert "ERR250683-tiny_bowtie2_aligned.sam" in sams[0].sam.path
    assert "ERR250683-tiny_bowtie2_report_aligned.txt" in sams[0].report.path

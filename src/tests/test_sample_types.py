from pathlib import Path
from flytekit import workflow
from flytekit.types.directory import FlyteDirectory

from config import test_assets
from tasks.utils import prepare_raw_samples
from tasks.sample_types import RawSample, FiltSample, SamFile


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


def test_sam_file_make_all():
    sams = SamFile.make_all(Path(test_assets["sam_dir"]))
    assert len(sams) == 2
    for i in sams:
        assert isinstance(i, SamFile)

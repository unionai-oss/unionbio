import os
import filecmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from tasks.sample_types import RawSample, FiltSample
from tasks.fastp import pyfastp
from tasks.fastqc import fastqc
from config import test_assets


def test_fastqc():
    qc_samp = fastqc(seq_dir=test_assets["seq_dir"])
    assert isinstance(qc_samp, FlyteDirectory)
    assert all(i in os.listdir(test_assets["fastqc_dir"]) for i in os.listdir(qc_samp.path))


def test_fastp():
    raw_samp = RawSample.make_all(Path(test_assets["seq_dir"]))[0]
    filt_samp = pyfastp(rs=raw_samp)
    assert isinstance(filt_samp, FiltSample)
    assert filecmp.cmp(
        Path(filt_samp.filt_r1.path),
        Path(test_assets["filt_dir"]).joinpath("ERR250683-tiny_1.filt.fastq.gz"),
    )

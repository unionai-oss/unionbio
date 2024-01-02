import os
import filecmp
from pathlib import Path
from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.utils import prepare_raw_samples
from tasks.sample_types import RawSample, FiltSample, SamFile
from tasks.fastp import pyfastp
from tasks.fastqc import fastqc
from config import test_assets, logger
from tests.utils import dir_contents_match


def test_fastqc():
    qc_samp = fastqc(seq_dir=test_assets["seq_dir"])
    print(os.listdir(qc_samp.path))
    assert isinstance(qc_samp, FlyteDirectory)
    assert filecmp.cmp(
        Path(qc_samp.path).joinpath("ERR250683-tiny_1_fastqc.html"),
        Path(test_assets["fastqc_dir"]).joinpath("ERR250683-tiny_1_fastqc.html"),
    )


def test_fastp():
    raw_samp = RawSample.make_all(Path(test_assets["seq_dir"]))[0]
    filt_samp = pyfastp(rs=raw_samp)
    assert isinstance(filt_samp, FiltSample)
    assert filecmp.cmp(
        Path(filt_samp.filt_r1.path),
        Path(test_assets["filt_dir"]).joinpath("ERR250683-tiny_1.filt.fastq.gz"),
    )

def test_multiqc():
    ...

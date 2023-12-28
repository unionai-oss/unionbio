from pathlib import Path
from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.utils import prepare_raw_samples
from tasks.sample_types import RawSample, FiltSample, SamFile
from tasks.fastp import pyfastp
from tasks.fastqc import fastqc
from config import test_assets, logger

def test_fastqc():
    qc_samp = fastqc(seq_dir=test_assets['local_tiny_seq_dir'])
    logger.debug(f'fastqc returned {qc_samp}')
    assert isinstance(qc_samp, FlyteDirectory)

def test_fastp():
    raw_samp = RawSample.make_all(Path(test_assets['local_seq_dir']))[0]
    filt_samp = pyfastp(rs=raw_samp)
    assert isinstance(filt_samp, FiltSample)

def test_multiqc():
    ...

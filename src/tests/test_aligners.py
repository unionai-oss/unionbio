from pathlib import Path
from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.utils import prepare_raw_samples
from tasks.hisat2 import hisat2_align_paired_reads, hisat2_index
from tasks.bowtie2 import bowtie2_align_samples, bowtie2_index
from config import test_assets, logger

def test_hisat2_index():
    idx_dir = hisat2_index(ref=Path(test_assets['local_tiny_ref']))
    assert isinstance(idx_dir, FlyteDirectory)

# def test_bowtie2_index():
#     ...
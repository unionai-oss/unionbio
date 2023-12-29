import os
from pathlib import Path
from filecmp import cmp
from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.utils import prepare_raw_samples
from tasks.hisat2 import hisat2_align_paired_reads, hisat2_index
from tasks.bowtie2 import bowtie2_align_samples, bowtie2_index
from config import test_assets, logger


def test_hisat2_index():
    idx_dir = hisat2_index(ref=Path(test_assets["ref_dir"]))
    assert isinstance(idx_dir, FlyteDirectory)
    assert all(
        x in os.listdir(test_assets["idx_dir"]) for x in os.listdir(idx_dir.path)
    )
    assert all(
        cmp(*i)
        for i in [
            (Path(test_assets["idx_dir"]).joinpath(x), Path(idx_dir.path).joinpath(x))
            for x in os.listdir(idx_dir.path)
        ]
    )

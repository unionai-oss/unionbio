import os
from pathlib import Path
from filecmp import cmp
from flytekit.types.directory import FlyteDirectory
from tasks.hisat2 import hisat2_index, hisat2_align_paired_reads
from tasks.bowtie2 import bowtie2_index, bowtie2_align_paired_reads
from datatypes.alignment import Alignment
from datatypes.reads import Reads
from config import test_assets


def test_hisat2_index():
    idx_dir = hisat2_index(ref=Path(test_assets["ref_path"]))
    assert isinstance(idx_dir, FlyteDirectory)
    assert all(
        x in os.listdir(test_assets["hs2_idx_dir"]) for x in os.listdir(idx_dir.path)
    )
    assert all(
        cmp(*i)
        for i in [
            (
                Path(test_assets["hs2_idx_dir"]).joinpath(x),
                Path(idx_dir.path).joinpath(x),
            )
            for x in os.listdir(idx_dir.path)
        ]
    )


def test_hisat2_align():
    idx_dir = FlyteDirectory(test_assets["hs2_idx_dir"])
    filt_samples = Reads.make_all(Path(test_assets["filt_dir"]))
    sam = hisat2_align_paired_reads(idx=idx_dir, fs=filt_samples[0])
    assert isinstance(sam, Alignment)
    assert all(
        x in os.listdir(test_assets["hs2_sam_dir"])
        for x in [Path(i).name for i in [sam.sam.path, sam.alignment_report.path]]
    )


def test_bowtie2_index():
    idx_dir = bowtie2_index(ref=Path(test_assets["ref_path"]))
    assert isinstance(idx_dir, FlyteDirectory)
    print(os.listdir(idx_dir.path))
    assert all(
        x in os.listdir(test_assets["bt2_idx_dir"]) for x in os.listdir(idx_dir.path)
    )
    assert all(
        cmp(*i)
        for i in [
            (
                Path(test_assets["bt2_idx_dir"]).joinpath(x),
                Path(idx_dir.path).joinpath(x),
            )
            for x in os.listdir(idx_dir.path)
        ]
    )


def test_bowtie2_align():
    idx_dir = FlyteDirectory(test_assets["bt2_idx_dir"])
    filt_samples = Reads.make_all(Path(test_assets["filt_dir"]))
    sam = bowtie2_align_paired_reads(idx=idx_dir, fs=filt_samples[0])
    assert isinstance(sam, Alignment)
    assert all(
        x in os.listdir(test_assets["bt2_sam_dir"])
        for x in [Path(i).name for i in [sam.sam.path, sam.alignment_report.path]]
    )

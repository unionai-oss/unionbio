import os
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from tasks.sample_types import RawSample, FiltSample, Alignment
from tasks.fastp import pyfastp
from tasks.fastqc import fastqc
from tasks.mark_dups import mark_dups
from tasks.sort_sam import sort_sam
from config import test_assets


def test_fastqc():
    qc_samp = fastqc(seq_dir=test_assets["seq_dir"])
    assert isinstance(qc_samp, FlyteDirectory)
    assert all(
        i in os.listdir(test_assets["fastqc_dir"]) for i in os.listdir(qc_samp.path)
    )


def test_fastp():
    raw_samp = RawSample.make_all(Path(test_assets["seq_dir"]))[0]
    filt_samp = pyfastp(rs=raw_samp)
    assert isinstance(filt_samp, FiltSample)
    assert cmp(
        Path(filt_samp.filt_r1.path),
        Path(test_assets["filt_dir"]).joinpath("ERR250683-tiny_1.filt.fastq.gz"),
    )


def test_sort_sam():
    alignment = Alignment.make_all(Path(test_assets["bt2_sam_dir"]))[0]
    alignment.sorted = True
    fname = alignment.get_alignment_fname()
    sorted_alignment = sort_sam(out_fname=fname, sam=alignment.sam)
    assert isinstance(sorted_alignment, FlyteFile)
    assert cmp(
        Path(sorted_alignment.path),
        Path(test_assets["sort_dir"]).joinpath(
            "ERR250683-tiny_bowtie2_sorted_aligned.sam"
        ),
    )


def test_mark_dups():
    alignment = Alignment.make_all(Path(test_assets["sort_dir"]))[0]
    alignment.deduped = True
    print(alignment)
    deduped, metrics = mark_dups(
        oafn=alignment.get_alignment_fname(),
        omfn=alignment.get_metrics_fname(),
        al=alignment.sam,
    )
    assert isinstance(deduped, FlyteFile)
    assert all(
        cmp(*i)
        for i in [
            (Path(test_assets["dedup_dir"]).joinpath(x), Path("/tmp/dedup").joinpath(x))
            for x in [deduped.path, metrics.path]
        ]
    )

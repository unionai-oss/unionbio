import os
import shutil
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.reference import Reference
from unionbio.tasks.fastp import pyfastp
from unionbio.tasks.fastqc import fastqc
from unionbio.tasks.bwa import bwa_index
from unionbio.tasks.mark_dups import mark_dups
from unionbio.tasks.sort_sam import sort_sam
from tests.config import test_assets
from tests.utils import dir_contents_match


def test_fastqc():
    qc_samp = fastqc(seq_dir=test_assets["raw_seq_dir"])
    assert isinstance(qc_samp, FlyteDirectory)
    assert all(
        i in os.listdir(test_assets["fastqc_dir"]) for i in os.listdir(qc_samp.path)
    )


def test_fastp():
    raw_samp = Reads.make_all(Path(test_assets["raw_seq_dir"]))[0]
    filt_samp = pyfastp(rs=raw_samp)
    assert isinstance(filt_samp, Reads)
    r1a = Path(filt_samp.read1.path)
    r1e = Path(test_assets["filt_seq_dir"]).joinpath("ERR250683-tiny_1.filt.fastq.gz")
    assert cmp(r1a, r1e)


def test_bwa_index(tmp_path):
    shutil.copy(test_assets["ref_path"], tmp_path)
    ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=tmp_path))
    indexed_ref = bwa_index(ref_obj=ref)
    assert isinstance(indexed_ref, Reference)
    assert f'{test_assets["ref_fn"]}.fai' in os.listdir(tmp_path)
    assert cmp(test_assets["ref_idx_path"], Path(tmp_path).joinpath(f'{test_assets["ref_fn"]}.fai'))
    assert dir_contents_match(test_assets["bwa_idx_dir"], tmp_path)


def test_sort_sam():
    alignment = Alignment.make_all(Path(test_assets["bt2_sam_dir"]))[0]
    alignment.sorted = True
    fname = alignment.get_alignment_fname()
    sorted_alignment = sort_sam(out_fname=fname, sam=alignment.alignment)
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
        al=alignment.alignment,
    )
    assert isinstance(deduped, FlyteFile)
    assert all(
        cmp(*i)
        for i in [
            (Path(test_assets["dedup_dir"]).joinpath(x), Path("/tmp/dedup").joinpath(x))
            for x in [deduped.path, metrics.path]
        ]
    )

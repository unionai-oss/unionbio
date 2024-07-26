import os
import shutil
from pathlib import Path
from filecmp import cmp
from flytekit.types.directory import FlyteDirectory
from unionbio.tasks.hisat2 import hisat2_index, hisat2_align_paired_reads
from unionbio.tasks.bowtie2 import bowtie2_index, bowtie2_align_paired_reads
from unionbio.tasks.bwa import bwa_index, bwa_align
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.reference import Reference
from tests.config import test_assets
from tests.utils import dir_conts_match, copy_dir_conts


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
    filt_samples = Reads.make_all(Path(test_assets["filt_seq_dir"]))
    al = hisat2_align_paired_reads(idx=idx_dir, fs=filt_samples[0])
    assert isinstance(al, Alignment)
    assert all(
        x in os.listdir(test_assets["hs2_sam_dir"])
        for x in [Path(i).name for i in [al.alignment.path, al.alignment_report.path]]
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
    filt_samples = Reads.make_all(Path(test_assets["filt_seq_dir"]))
    al = bowtie2_align_paired_reads(idx=idx_dir, fs=filt_samples[0])
    assert isinstance(al, Alignment)
    assert all(
        x in os.listdir(test_assets["bt2_sam_dir"])
        for x in [Path(i).name for i in [al.alignment.path, al.alignment_report.path]]
    )


def test_bwa_index(tmp_path):
    copy_dir_conts(test_assets["ref_dir"], tmp_path)
    ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=tmp_path))
    indexed_ref = bwa_index(ref_obj=ref)
    assert isinstance(indexed_ref, Reference)
    assert f'{test_assets["ref_fn"]}.fai' in os.listdir(tmp_path)
    assert cmp(test_assets["ref_idx_path"], Path(tmp_path).joinpath(f'{test_assets["ref_fn"]}.fai'))
    assert dir_conts_match(test_assets["bwa_idx_dir"], tmp_path)


def test_bwa_align(tmp_path):
    copy_dir_conts(test_assets["filt_seq_dir"], tmp_path)
    reads = Reads.make_all(tmp_path)[0]
    shutil.copy(test_assets["ref_path"], tmp_path)
    copy_dir_conts(test_assets["bwa_idx_dir"], tmp_path)
    ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=tmp_path), test_assets["ref_fn"], "bwa")
    print("IN TEST:")
    print(ref.ref_dir.path)
    print(reads.read1.path)
    alignment = bwa_align(ref=ref, reads=reads)
    # assert isinstance(alignment, Alignment)
    # assert cmp(
    #     Path(alignment.get_alignment_fname()),
    #     Path(test_assets["bwa_sam_dir"]).joinpath("ERR250683-tiny_bwa_aligned.sam"),
    # )
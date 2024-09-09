import os
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from unionbio.tasks.fastp import pyfastp
from unionbio.tasks.fastqc import fastqc
from unionbio.tasks.mark_dups import mark_dups
from unionbio.tasks.sort_sam import sort_sam
from unionbio.tasks.base_recal import base_recalibrator
from tests.config import test_assets
from tests.utils import copy_dir_conts
from unionbio.types import Alignment, Reads, Reference, VCF


def test_fastqc(tmp_path):
    copy_dir_conts(test_assets["raw_seq_dir"], tmp_path)
    reads = Reads.make_all(tmp_path)
    qc_dir = fastqc(reads=reads)
    assert isinstance(qc_dir, FlyteDirectory)
    assert all(
        i in os.listdir(test_assets["fastqc_dir"]) for i in os.listdir(qc_dir.path)
    )


def test_fastp():
    raw_samp = Reads.make_all(test_assets["raw_seq_dir"])[0]
    filt_samp = pyfastp(rs=raw_samp)
    assert isinstance(filt_samp, Reads)
    r1a = Path(filt_samp.read1.path)
    r1e = test_assets["filt_seq_dir"].joinpath("ERR250683-tiny_1.filt.fastq.gz")
    assert cmp(r1a, r1e)


def test_sort_sam():
    alignment = Alignment.make_all(test_assets["bt2_sam_dir"])[0]
    sorted_alignment = sort_sam(al=alignment)
    assert isinstance(sorted_alignment, Alignment)


def test_mark_dups():
    alignment = Alignment.make_all(test_assets["sort_dir"])[0]
    dd_al = mark_dups(al=alignment)
    assert isinstance(dd_al, Alignment)
    assert dd_al.deduped
    assert Path(dd_al.alignment.path).exists()


def test_base_recal(tmp_path):
    copy_dir_conts(test_assets["dedup_dir"], tmp_path)
    copy_dir_conts(test_assets["ref_dir"], tmp_path)
    copy_dir_conts(test_assets["sites_dir"], tmp_path)
    alignment = Alignment.make_all(tmp_path, include=["SRR812824-chr21*"])[0]
    ref = Reference("GRCh38_chr21.fasta", FlyteDirectory(path=tmp_path))
    sites = VCF.make_all(tmp_path, include=["*gold_standard_chr21*"])[0]
    recal_out = base_recalibrator(al=alignment, ref=ref, sites=sites)
    assert isinstance(recal_out, Alignment)
    assert recal_out.recalibrated

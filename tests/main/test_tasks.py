import os
import shutil
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.tasks.fastp import pyfastp
from unionbio.tasks.fastqc import fastqc
from unionbio.tasks.bwa import bwa_index, bwa_align
from unionbio.tasks.mark_dups import mark_dups
from unionbio.tasks.sort_sam import sort_sam
from unionbio.tasks.base_recal import recalibrate_bases
from tests.config import test_assets
from tests.utils import dir_conts_match, copy_dir_conts, comp_files


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


def test_sort_sam():
    alignment = Alignment.make_all(Path(test_assets["bt2_sam_dir"]))[0]
    sorted_alignment = sort_sam(al=alignment)
    assert isinstance(sorted_alignment, Alignment)


def test_mark_dups():
    alignment = Alignment.make_all(Path(test_assets["sort_dir"]))[0]
    dd_al = mark_dups(al=alignment)
    assert isinstance(dd_al, Alignment)
    assert dd_al.deduped
    assert Path(dd_al.alignment.path).exists()

def test_base_recal(tmp_path):
    copy_dir_conts(test_assets["recal_in"], tmp_path)
    copy_dir_conts(Path(test_assets["ref_dir"]).joinpath("chr21"), tmp_path)
    copy_dir_conts(Path(test_assets["sites_dir"]).joinpath("chr21"), tmp_path)
    alignment = Alignment.make_all(tmp_path)[0]
    ref = Reference("GRCh38_chr21.fasta", FlyteDirectory(path=tmp_path))
    sites = VCF.make_all(tmp_path)[0]
    recal_out = recalibrate_bases(al=alignment, ref=ref, sites=sites)
    assert isinstance(recal_out, Alignment)
    assert recal_out.recalibrated
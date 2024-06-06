import os
import string
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.variants import VCF
from unionbio.tasks.utils import fetch_file, intersect_vcfs
from unionbio.tasks.helpers import gunzip_file
from tests.config import test_assets


def test_fetch_http_file(tmp_path):
    # Test that fetch_file downloads a file
    url = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
    outpath = fetch_file(url, tmp_path)
    print(outpath)
    assert outpath.exists()
    assert open(outpath, "rb").read(), "File {outpath} is empty"
    assert outpath.name == "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"


def test_fetch_ftp_file(tmp_path):
    # Test that fetch_file downloads a file
    url = "ftp://ftp.ncbi.nlm.nih.gov/tech-reports/tech-report.txt"
    outpath = fetch_file(url, tmp_path)
    assert outpath.exists()
    assert open(outpath, "rb").read(), "File {outpath} is empty"
    assert outpath.name == "tech-report.txt"


def test_intersect_vcfs():
    vcf1 = VCF.make_all(Path(test_assets["vcf_dir"]))[0]
    vcf2 = VCF.make_all(Path(test_assets["vcf_dir"]))[0]
    out = intersect_vcfs(vcf1=vcf1, vcf2=vcf2)
    print(out)
    assert isinstance(out, VCF)


def test_gunzip():
    gzfile = test_assets["sites_path"]
    unzipped = gunzip_file(Path(gzfile))
    assert unzipped.exists()
    assert all([c in string.printable for c in open(unzipped).readline().strip()])

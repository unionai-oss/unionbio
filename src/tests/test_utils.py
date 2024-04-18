import os
import string
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from datatypes.alignment import Alignment
from datatypes.reads import Reads
from tasks.utils import fetch_file, intersect_vcfs
from tasks.helpers import gunzip_file
from config import test_assets


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


def test_intersect_vcfs(tmp_path):
    # Test that intersect_vcfs returns a FlyteFile
    vcf1 = tmp_path / "vcf1.vcf"
    vcf2 = tmp_path / "vcf2.vcf"
    vcf1.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    vcf2.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    vcf1_obj = FlyteFile(vcf1)
    vcf2_obj = FlyteFile(vcf2)
    out = intersect_vcfs(vcf1=vcf1_obj, vcf2=vcf2_obj)
    assert isinstance(out, FlyteFile)

def test_gunzip():
    gzfile = test_assets['sites_path']
    unzipped = gunzip_file(Path(gzfile))
    assert unzipped.exists()
    assert all([c in string.printable for c in open(unzipped).readline().strip()])

import string
from pathlib import Path
from filecmp import cmp
from unionbio.datatypes.variants import VCF
from unionbio.tasks.utils import intersect_vcfs
from unionbio.tasks.helpers import gunzip_file, fetch_file
from tests.config import test_assets


def test_fetch_http_file(tmp_path):
    # Test that fetch_file downloads a file
    url = "https://github.com/unionai-oss/unionbio/raw/main/tests/assets/sequences/raw/ERR250683-tiny_1.fastq.gz"
    outpath = fetch_file(url, tmp_path)
    assert outpath.exists()
    assert open(outpath, "rb").read(), "File {outpath} is empty"
    assert outpath.name == "ERR250683-tiny_1.fastq.gz"
    assert cmp(outpath, Path(test_assets["raw_read"]))

def test_fetch_ftp_file(tmp_path):
    # Test that fetch_file downloads a file
    url = "ftp://ftp.ncbi.nlm.nih.gov/tech-reports/tech-report.txt"
    outpath = fetch_file(url, tmp_path)
    assert outpath.exists()
    assert open(outpath, "rb").read(), "File {outpath} is empty"
    assert outpath.name == "tech-report.txt"


def test_intersect_vcfs():
    vcfs = VCF.make_all(Path(test_assets["vcf_dir"]))
    # vcf2 = VCF.make_all(Path(test_assets["vcf_dir"]))[0]
    out = intersect_vcfs(vcf1=vcfs[0], vcf2=vcfs[1])
    print(out)
    assert isinstance(out, VCF)


def test_gunzip():
    gzfile = test_assets["sites_path"]
    unzipped = gunzip_file(Path(gzfile))
    assert unzipped.exists()
    assert all([c in string.printable for c in open(unzipped).readline().strip()])

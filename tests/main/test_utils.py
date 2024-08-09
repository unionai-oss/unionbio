import string
import shutil
from pathlib import Path
from filecmp import cmp
from unionbio.datatypes.variants import VCF
from unionbio.datatypes.alignment import Alignment
from unionbio.tasks.utils import intersect_vcfs, reformat_alignments
from unionbio.tasks.helpers import gunzip_file, fetch_file
from tests.config import test_assets
from tests.utils import copy_dir_conts


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

def test_intersect_vcfs(tmp_path):
    copy_dir_conts(test_assets["vcf_dir"], tmp_path)
    vcfs = VCF.make_all(tmp_path)
    out = intersect_vcfs(vcf1=vcfs[0], vcf2=vcfs[1])
    out.aggregate(target=tmp_path)
    assert isinstance(out, VCF)

def test_gunzip():
    gzfile = test_assets["sites_path"]
    unzipped = gunzip_file(Path(gzfile))
    assert unzipped.exists()
    assert all([c in string.printable for c in open(unzipped).readline().strip()])

def test_alignment_reformat(tmp_path):
    shutil.copy2(Path(test_assets["al_dir"]).joinpath("recal/SRR812824-chr21_bowtie2_sorted_deduped_recal.sam"), tmp_path)
    al_in = Alignment.make_all(tmp_path)
    exp_out = Alignment.make_all(Path(test_assets["bam_dir"]))[0]
    al_out = reformat_alignments(al_in, to_format="bam")[0]
    ep = Path(exp_out.alignment.path)
    ap = Path(al_out.alignment.path)
    assert isinstance(al_out, Alignment)
    assert ap.exists()
    assert ep.name == ap.name


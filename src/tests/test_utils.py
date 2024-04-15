import os
from filecmp import cmp
from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from datatypes.alignment import Alignment
from datatypes.reads import Reads
from tasks.utils import fetch_file, fetch_remote_reads, fetch_remote_reference, fetch_remote_sites, intersect_vcfs

def test_fetch_file(tmp_path):
    # Test that fetch_file downloads a file
    url = 'https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
    outpath = fetch_file(url, tmp_path)
    print(outpath)
    assert outpath.exists()
    assert open(outpath, 'rb').read(), 'File {outpath} is empty'
    assert outpath.name == 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'

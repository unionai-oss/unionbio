from pathlib import Path
from flytekit.types.directory import FlyteDirectory
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.tasks.haplotype_caller import haplotype_caller
from tests.utils import copy_dir_conts
from tests.config import test_assets

def test_call_variants(tmp_path):
    copy_dir_conts(Path(test_assets["ref_dir"]).joinpath("chr21"), tmp_path)
    copy_dir_conts(test_assets["bam_dir"], tmp_path)
    ref = Reference("GRCh38_chr21.fasta", FlyteDirectory(path=tmp_path))
    al = Alignment.make_all(tmp_path)[0]
    vcf = haplotype_caller(ref=ref, al=al)
    assert isinstance(vcf, VCF)
from flytekit.types.directory import FlyteDirectory
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.tasks.haplotype_caller import hc_call_variants
from tests.utils import copy_dir_conts
from tests.config import test_assets

def test_call_variants():#(tmp_path):
    tmp_path = "/tmp/hc_call_variants"
    copy_dir_conts(test_assets["ref_dir"], tmp_path)
    copy_dir_conts(test_assets["sort_dir"], tmp_path)
    # ref = Reference(test_assets["ref_fn"], FlyteDirectory(path=tmp_path))
    # al = Alignment.make_all(tmp_path)[0]
    # vcf = hc_call_variants(ref=ref, al=al)
    # assert isinstance(vcf, VCF)
from tasks.parabricks import pb_fq2bam, pb_deepvar, pb_haplocall
from datatypes.reference import Reference
from datatypes.reads import Reads
from datatypes.alignment import Alignment
from datatypes.variants import VCF
from config import test_assets

ref_obj = Reference(ref_name=test_assets["ref_fn"], ref_dir=test_assets["ref_dir"])
sites_obj = VCF(vcf=test_assets["sites_path"], vcf_idx=test_assets["sites_idx_path"])

def test_pb_fq2bam():
    reads_obj = Reads.make_all(test_assets["seq_dir"])
    alignment = pb_fq2bam(reads=reads_obj, sites=sites_obj, ref=ref_obj)
    assert isinstance(alignment, Alignment)

def test_pb_deepvar():
    al_obj = Alignment.make_all(test_assets["pb_fq2bam_dir"])
    deepvar_vcf = pb_deepvar(al=al_obj, ref=ref_obj)
    assert isinstance(deepvar_vcf, VCF)

def test_pb_haplocall():
    al_obj = Alignment.make_all(test_assets["pb_haplocall_dir"])
    haplocall_vcf = pb_haplocall(al=al_obj, ref=ref_obj)
    assert isinstance(haplocall_vcf, VCF)

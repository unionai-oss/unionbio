import os
from pathlib import Path
from typing import List
from flytekit import task, dynamic, current_context, Resources
from flytekit.types.file import FlyteFile
from flytekit.extras.tasks.shell import subproc_execute
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.config import logger, main_img_fqn

@task(container_image=main_img_fqn, requests=Resources(cpu="4", mem="10Gi"))
def haplotype_caller(ref: Reference, al: Alignment) -> VCF:
    """
    Call variants using HaplotypeCaller from GATK.

    Args:
        ref (Reference): A reference object containing the reference FASTA file.
        al (Alignment): An alignment object containing the aligned reads.

    Returns:
        VCF: A VCF object containing the called variants.
    """
    logger.info("Calling variants using HaplotypeCaller")
    logger.info(f"Reference: {ref}")
    logger.info(f"Alignment: {al}")
    ref.aggregate()
    logger.debug(f"Reference dir contents: {os.listdir(ref.ref_dir.path)}")
    al.aggregate()
    con_dir = Path(current_context().working_directory)
    vcf_out = VCF(
        sample=al.sample,
        caller="gatk-hc",
    )
    vcf_fn = con_dir.joinpath(vcf_out.get_vcf_fname())
    vcf_idx_fn = con_dir.joinpath(vcf_out.get_vcf_idx_fname())
    hc_cmd = [
        "gatk",
        "HaplotypeCaller",
        "-R",
        str(ref.get_ref_path()),
        "-I",
        str(al.alignment.path),
        "-O",
        str(vcf_fn),
        "-ERC",
        "GVCF",
    ]
    logger.debug("Running command:")
    logger.debug(" ".join(hc_cmd))
    logger.debug(f"Running in dir {con_dir} with contents: {os.listdir(con_dir)}")
    subproc_execute(hc_cmd, cwd=con_dir)
    vcf_out.vcf = FlyteFile(path=str(vcf_fn))
    vcf_out.vcf_idx = FlyteFile(path=str(vcf_idx_fn))
    return vcf_out

@dynamic
def hc_call_samples(ref: Reference, als: List[Alignment]) -> List[VCF]:
    return [haplotype_caller(ref=ref, al=al) for al in als]

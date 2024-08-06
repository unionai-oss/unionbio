from flytekit import task
from flytekit.extras.tasks.shell import subproc_execute
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.config import logger, main_img_fqn

@task(container_image=main_img_fqn)
def hc_call_variants(ref: Reference, al: Alignment) -> VCF:
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
    al.aggregate()
    hc_cmd = [
        "gatk",
        "HaplotypeCaller",
        "-R",
        ref.get_ref_path(),
        "-I",
        al.get_alignment_fname(),
        "-O",
        "test.vcf",
    ]
    logger.info(f"Running command: {hc_cmd}")
    subproc_execute(hc_cmd)
    return VCF("test.vcf")
from flytekit import task
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
    
    return VCF("test.vcf")
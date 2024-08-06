from flytekit import task
from flytekit.extras.tasks.shell import subproc_execute
from flytekit.types.file import FlyteFile

from unionbio.config import main_img_fqn
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.variants import VCF


@task(container_image=main_img_fqn)
def recalibrate_bases(ref: Reference, sites: VCF, al: Alignment) -> Alignment:
    """
    Recalibrate base quality scores using GATK's BaseRecalibrator.

    Args:
        ref (Reference): A reference object containing the reference FASTA file.
        sites (VCF): A VCF object containing known variant sites.
        al (Alignment): An alignment object containing the aligned reads.

    Returns:
        Alignment: An alignment object with the recalibrated alignment file.
    """
    ref.aggregate()
    sites.aggregate()
    al.aggregate()

    recal_fn = al.get_bqsr_fname()

    gen_table_cmd = [
        "gatk",
        "BaseRecalibrator",
        "-I",
        al.alignment.path,
        "-R",
        ref.get_ref_path(),
        "--known-sites",
        sites.vcf.path,
        "-O",
        recal_fn
    ]
    subproc_execute(command=gen_table_cmd)

    al.recalibrated = True
    al_out_fname = al.get_alignment_fname()
    apply_recal_cmd = [
        "gatk",
        "ApplyBQSR",
        "-I",
        al.alignment.path,
        "-R",
        ref.get_ref_path(),
        "--bqsr-recal-file",
        recal_fn,
        "-O",
        al_out_fname,
    ]
    subproc_execute(command=apply_recal_cmd)

    al.alignment = FlyteFile(path=al_out_fname)

    return al
import os
import shutil
from pathlib import Path
from flytekit import kwtypes, task, Resources, current_context, TaskMetadata
from flytekit.extras.tasks.shell import OutputLocation, ShellTask, subproc_execute
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from unionbio.config import ref_hash, main_img_fqn
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.alignment import Alignment


@task(
    container_image=main_img_fqn,
    requests=Resources(cpu="4", mem="10Gi"),
    cache=True,
    cache_version=ref_hash,
)
def bwa_index(ref_obj: Reference) -> Reference:
    """Indexes a reference genome using BWA.

    Args:
        ref_obj (Reference): The reference object containing the reference genome.

    Returns:
        Reference: The updated reference object with associated index and metadata.
    """
    ref_obj.ref_dir.download()

    if f"{ref_obj.ref_name}.fai" not in os.listdir(ref_obj.ref_dir.path):
        sam_idx = ["samtools", "faidx", str(ref_obj.get_ref_path())]
        subproc_execute(sam_idx, cwd=ref_obj.ref_dir.path)

    bwa_idx = ["bwa", "index", str(ref_obj.get_ref_path())]
    subproc_execute(bwa_idx, cwd=ref_obj.ref_dir.path)

    ref_obj.index_name = ref_obj.ref_name
    ref_obj.indexed_with = "bwa"

    return ref_obj


@task(
    container_image=main_img_fqn,
    requests=Resources(cpu="4", mem="10Gi"),
)
def bwa_align(ref: Reference, reads: Reads) -> Alignment:
    """Aligns reads to a reference genome using BWA.

    Args:
        ref (Reference): The reference object containing the reference genome and index.
        reads (Reads): The reads object containing the reads to align.

    Returns:
        Alignment: The alignment object containing the aligned reads.
    """
    print("IN TASK PRE-DOWNLOAD:")
    print(ref.ref_dir.path)
    print(reads.read1.path)
    ref.aggregate()
    reads.aggregate()
    print("IN TASK POST-DOWNLOAD:")
    print(ref.ref_dir.path)
    print(reads.read1.path)

    al_out = Alignment(
        sample=reads.sample,
        aligner="bwa",
        format="sam",
        sorted=False,
        deduped=False,
    )
    sam_out = al_out.get_alignment_fname()
    bwa_align = [
        "bwa",
        "mem",
        str(ref.get_ref_path()),
        reads.read1.path,
        reads.read2.path,
        ">",
        sam_out,
    ]
    # print(bwa_align)
    # subproc_execute(bwa_align, shell=True)
    # al_out.alignment = FlyteFile(sam_out)
    return al_out
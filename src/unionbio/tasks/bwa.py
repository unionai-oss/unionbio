import os
import shutil
import subprocess
from pathlib import Path
from typing import List
from flytekit import kwtypes, task, Resources, current_context, TaskMetadata, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask, subproc_execute
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from unionbio.config import remote_ref, main_img_fqn, logger
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.alignment import Alignment


@task(
    container_image=main_img_fqn,
    requests=Resources(cpu="4", mem="10Gi"),
    cache=True,
    cache_version=remote_ref,
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
    ref.aggregate()
    reads.aggregate()

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
        str(reads.read1.path),
        str(reads.read2.path),
        ">",
        sam_out,
    ]
    cmd_str = " ".join(bwa_align)
    logger.info(f"Running BWA alignment with: {cmd_str}")
    subproc_execute(cmd_str, shell=True)
    sp = Path(sam_out)
    logger.debug(f"Alignment exists ({sp.exists()}) at {sp.resolve()}")
    al_out.alignment = FlyteFile(path=sam_out)
    return al_out

@dynamic(container_image=main_img_fqn)
def bwa_align_samples(idx: Reference, samples: List[Reads]) -> List[Alignment]:
    """
    Process samples through BWA.

    This function takes a Reference object containing a bwa index and a list of
    Reads objects containing filtered sample data. It performs paired-end alignment
    using bwa. It then returns a list of Alignment objects representing the alignment results.

    Args:
        idx (FlyteDirectory): The FlyteDirectory object representing the bowtie2 index.
        samples (List[Reads]): A list of Reads objects containing sample data
            to be processed.

    Returns:
        [List[Alignment]: A list of alignment results (Alignment objects).
    """
    return [bwa_align(ref=idx, reads=sample) for sample in samples]
from pathlib import Path
from typing import List
from flytekit import kwtypes, task, Resources, current_context, TaskMetadata, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask, subproc_execute
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from unionbio.config import ref_hash, main_img_fqn, logger
from unionbio.datatypes.alignment import Alignment
from unionbio.datatypes.reads import Reads


"""
Generate Bowtie2 index files from a reference genome.

Args:
    ref (FlyteFile): A FlyteFile object representing the input reference file.

Returns:
    FlyteDirectory: A FlyteDirectory object containing the index files.
"""
bowtie2_index = ShellTask(
    name="bowtie2-index",
    debug=True,
    requests=Resources(cpu="4", mem="10Gi"),
    metadata=TaskMetadata(retries=3, cache=True, cache_version=ref_hash),
    container_image=main_img_fqn,
    script="""
    mkdir {outputs.idx}
    bowtie2-build {inputs.ref} {outputs.idx}/bt2_idx
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[
        OutputLocation(var="idx", var_type=FlyteDirectory, location="/tmp/bt2_idx")
    ],
)


@task(
    container_image=main_img_fqn,
    requests=Resources(cpu="4", mem="10Gi"),
)
def bowtie2_align_paired_reads(idx: FlyteDirectory, fs: Reads) -> Alignment:
    """
    Perform paired-end alignment using Bowtie 2 on a filtered sample.

    This function takes a FlyteDirectory object representing the Bowtie 2 index and a
    FiltSample object containing filtered sample data. It performs paired-end alignment
    using Bowtie 2 and returns a Alignment object representing the resulting alignment.

    Args:
        idx (FlyteDirectory): A FlyteDirectory object representing the Bowtie 2 index.
        fs (Reads): A filtered sample Reads object containing filtered sample data to be aligned.

    Returns:
        Alignment: An Alignment object representing the alignment result.
    """
    idx.download()
    logger.debug(f"Index downloaded to {idx.path}")
    ldir = Path(current_context().working_directory)

    alignment = Alignment(fs.sample, "bowtie2", "sam")
    al = ldir.joinpath(alignment.get_alignment_fname())
    rep = ldir.joinpath(alignment.get_report_fname())
    logger.debug(f"Writing alignment to {al} and report to {rep}")

    cmd = [
        "bowtie2",
        "-x",
        f"{idx.path}/bt2_idx",
        "-1",
        fs.read1,
        "-2",
        fs.read2,
        "-S",
        al,
    ]
    logger.debug(f"Running command: {cmd}")

    result = subproc_execute(cmd)

    with open(rep, "w") as f:
        f.write(result.error)

    setattr(alignment, "alignment", FlyteFile(path=str(al)))
    setattr(alignment, "alignment_report", FlyteFile(path=str(rep)))
    setattr(alignment, "sorted", False)
    setattr(alignment, "deduped", False)

    return alignment


@dynamic(container_image=main_img_fqn)
def bowtie2_align_samples(idx: FlyteDirectory, samples: List[Reads]) -> List[Alignment]:
    """
    Process samples through bowtie2.

    This function takes a FlyteDirectory objects representing a bowtie index and a list of
    Reads objects containing filtered sample data. It performs paired-end alignment
    using bowtie2. It then returns a list of Alignment objects representing the alignment results.

    Args:
        bt2_idx (FlyteDirectory): The FlyteDirectory object representing the bowtie2 index.
        samples (List[Reads]): A list of Reads objects containing sample data
            to be processed.

    Returns:
        List[List[Alignment]]: A list of lists, where each inner list contains alignment
            results (Alignment objects) for a sample, with results from both aligners.
    """
    sams = []
    for sample in samples:
        sam = bowtie2_align_paired_reads(idx=idx, fs=sample)
        sams.append(sam)
    return sams

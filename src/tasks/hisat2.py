import gzip
import shutil
from pathlib import Path
from flytekit import kwtypes, task, Resources, current_context, TaskMetadata
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from config import ref_hash, base_image, logger
from datatypes.alignment import Alignment
from datatypes.reads import Reads
from tasks.utils import subproc_raise

"""
Generate Hisat2 index files from a reference genome.

Args:
    ref (FlyteFile): A FlyteFile object representing the input reference file.

Returns:
    FlyteDirectory: A FlyteDirectory object containing the index files.
"""
hisat2_index = ShellTask(
    name="hisat2-index",
    debug=True,
    metadata=TaskMetadata(retries=3, cache=True, cache_version=ref_hash),
    requests=Resources(cpu="4", mem="10Gi"),
    container_image=base_image,
    script="""
    mkdir {outputs.idx}
    hisat2-build {inputs.ref} {outputs.idx}/hs2_idx
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[
        OutputLocation(var="idx", var_type=FlyteDirectory, location="/tmp/hs2_idx")
    ],
)


@task(
    container_image=base_image,
    requests=Resources(cpu="4", mem="10Gi"),
)
def hisat2_align_paired_reads(idx: FlyteDirectory, fs: Reads) -> Alignment:
    """
    Perform paired-end alignment using Hisat 2 on a filtered sample.

    This function takes a FlyteDirectory object representing the Hisat 2 index and a
    Reads object containing filtered sample data. It performs paired-end alignment
    using Hisat 2 and returns a Alignment object representing the resulting alignment.

    Args:
        idx (FlyteDirectory): A FlyteDirectory object representing the Hisat 2 index.
        fs (Reads): A Reads object containing filtered sample data to be aligned.

    Returns:
        Alignment: An Alignment object representing the alignment result in SAM format.
    """
    idx.download()
    ldir = Path(current_context().working_directory)
    alignment = Alignment(fs.sample, "hisat2")
    sam = ldir.joinpath(alignment.get_alignment_fname())
    rep = ldir.joinpath(alignment.get_report_fname())
    logger.debug(f"Writing SAM to {sam} and report to {rep}")

    unc_r1 = ldir.joinpath(f"{fs.sample}_1.fq")
    unc_r2 = ldir.joinpath(f"{fs.sample}_2.fq")

    with gzip.open(fs.read1, "rb") as gz_file, open(unc_r1, "wb") as out_file:
        shutil.copyfileobj(gz_file, out_file)
    with gzip.open(fs.read2, "rb") as gz_file, open(unc_r2, "wb") as out_file:
        shutil.copyfileobj(gz_file, out_file)
    logger.debug(f"Uncompressed reads to {unc_r1} and {unc_r2}")

    cmd = [
        "hisat2",
        "-x",
        f"{idx.path}/hs2_idx",
        "-1",
        unc_r1,
        "-2",
        unc_r2,
        "-S",
        sam,
        "--summary-file",
        rep,
    ]
    logger.debug(f"Running command: {cmd}")

    stdout, stderr = subproc_raise(cmd)

    setattr(alignment, "sam", FlyteFile(path=str(sam)))
    setattr(alignment, "alignment_report", FlyteFile(path=str(rep)))
    setattr(alignment, "sorted", False)
    setattr(alignment, "deduped", False)

    return alignment

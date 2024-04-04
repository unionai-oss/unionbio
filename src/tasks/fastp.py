from pathlib import Path
from flytekit import task, Resources, current_context
from flytekit.types.file import FlyteFile

from config import base_image, logger, fastp_cpu
from datatypes.reads import Reads
from tasks.utils import subproc_raise


@task(
    requests=Resources(cpu=fastp_cpu, mem="2Gi"),
    container_image=base_image,
)
def pyfastp(rs: Reads) -> Reads:
    """
    Perform quality filtering and preprocessing using Fastp on a RawSample.

    This function takes a RawSample object containing raw sequencing data, performs quality
    filtering and preprocessing using the pyfastp tool, and returns a FiltSample object
    representing the filtered and processed data.

    Args:
        rs (RawSample): A RawSample object containing raw sequencing data to be processed.

    Returns:
        FiltSample: A FiltSample object representing the filtered and preprocessed data.
    """
    ldir = Path(current_context().working_directory)
    samp = Reads(rs.sample)
    samp.filtered = True
    o1, o2 = samp.get_read_fnames()
    rep = samp.get_report_fname()
    o1p = ldir.joinpath(o1)
    o2p = ldir.joinpath(o2)
    repp = ldir.joinpath(rep)
    logger.debug(f"Writing filtered reads to {o1p} and {o2p} and report to {repp}")

    cmd = [
        "fastp",
        "-i",
        rs.read1,
        "-I",
        rs.read2,
        "--thread",
        str(int(fastp_cpu) * 2),
        "-o",
        o1p,
        "-O",
        o2p,
        "-j",
        rep,
    ]
    logger.debug(f"Running command: {cmd}")

    subproc_raise(cmd)

    setattr(samp, "read1", FlyteFile(path=str(o1p)))
    setattr(samp, "read2", FlyteFile(path=str(o2p)))
    setattr(samp, "filt_report", FlyteFile(path=str(repp)))

    return samp

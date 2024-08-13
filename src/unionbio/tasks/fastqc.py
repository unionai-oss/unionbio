from pathlib import Path
from typing import List
from flytekit import kwtypes, TaskMetadata, task, current_context
from flytekit.extras.tasks.shell import OutputLocation, ShellTask, subproc_execute
from flytekit.types.directory import FlyteDirectory
from unionbio.config import main_img_fqn, logger
from unionbio.datatypes.reads import Reads

@task
def fastqc(reads: List[Reads]) -> FlyteDirectory:
    """
    Perform quality control using FastQC.

    This function takes a FlyteDirectory object containing raw sequencing data, 
    gathers QC metrics using FastQC, and returns a FlyteDirectory object that
    can be crawled with MultiQC to generate a report.

    Args:
        reads (List[Reads]): A list of Reads objects containing raw sequencing data.

    Returns:
        qc (FlyteDirectory): A directory containing fastqc report output.
    """
    pwd = Path(current_context().working_directory)
    indir = pwd.joinpath("fastqc_in")
    outdir = pwd.joinpath("fastqc_out")
    outdir.mkdir()
    for r in reads:
        r.aggregate(target=indir)
    
    fqc_cmd = [
        "fastqc",
        f"{indir}/*.fastq.gz",
        "--outdir",
        str(outdir),
    ]
    logger.info(f"Running FastQC with command: {fqc_cmd}")
    subproc_execute(fqc_cmd)
    
    return FlyteDirectory(path=outdir)

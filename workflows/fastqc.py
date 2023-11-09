from flytekit import kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory

from .config import base_image

"""
Perform quality control using FastQC.

This function takes a FlyteDirectory object containing raw sequencing data, 
gathers QC metrics using FastQC, and returns a FlyteDirectory object that
can be crawled with MultiQC to generate a report.

Args:
    seq_dir (FlyteDirectory): An S3 prefix containing raw sequencing data to be processed.

Returns:
    qc (FlyteDirectory): A directory containing fastqc report output.
"""
fastqc = ShellTask(
    name="fastqc",
    debug=True,
    script="""
    mkdir {outputs.qc}
    fastqc {inputs.seq_dir}/*.fastq.gz --outdir={outputs.qc}
    """,
    inputs=kwtypes(seq_dir=FlyteDirectory),
    output_locs=[
        OutputLocation(var="qc", var_type=FlyteDirectory, location="/root/qc")
    ],
    container_image=base_image,
)

from flytekit import kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory

from .config import base_image

# human in the loop after fastqc?
# surafce report in flytedeck
# or just have conditional for quality counts before proceeding
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

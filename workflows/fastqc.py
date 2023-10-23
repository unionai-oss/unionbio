from flytekit import kwtypes, task, workflow, ImageSpec
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

fastqc = ShellTask(
    name="fastqc",
    debug=True,
    script=
    """
    mkdir {outputs.qc}
    fastqc {inputs.seq_dir}/*.fastq.gz --outdir={outputs.qc}
    """,
    inputs=kwtypes(seq_dir=FlyteDirectory),
    output_locs=[OutputLocation(var="qc", var_type=FlyteDirectory, location='/root/qc')],
    container_image='localhost:30000/variant-discovery:latest'
)

@workflow
def fastqc_wf() -> FlyteDirectory:
    return fastqc(seq_dir='s3://my-s3-bucket/my-data/sequences')
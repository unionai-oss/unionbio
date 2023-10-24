from flytekit import kwtypes, task, workflow, ImageSpec, Resources
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

fastp = ShellTask(
    name="fastp",
    debug=True,
    requests=Resources(cpu="1", mem="2Gi"),
    script=
    """
    fastp -i {inputs.i1} -I {inputs.i2} -o {outputs.o1} -O {outputs.o2} -j {outputs.rep}
    """,
    inputs=kwtypes(i1=FlyteFile, i2=FlyteFile),
    output_locs=[
        OutputLocation(var="o1", var_type=FlyteFile, location='/root/out1.fq.gz'),
        OutputLocation(var="o2", var_type=FlyteFile, location='/root/out2.fq.gz'),
        OutputLocation(var="rep", var_type=FlyteFile, location='/root/fastp.json'),
        ],
    container_image='localhost:30000/variant-discovery:latest'
)

@workflow
def fastp_wf():
    fastp(
        i1='s3://my-s3-bucket/my-data/sequences/ERR250683_1.fastq.gz',
        i2='s3://my-s3-bucket/my-data/sequences/ERR250683_2.fastq.gz',
    )
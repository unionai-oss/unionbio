from flytekit import kwtypes, task, workflow, ImageSpec
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
# from workflows import config, logger

bowtie_image_spec = ImageSpec(
    base_image="refbase:latest",
    apt_packages=["bowtie2"]
)

check_base_image = ShellTask(
    name="check-base-image",
    debug=True,
    script=
    """
    ls -lah /root/workflows
    """,
    inputs=kwtypes(),
    output_locs=[],
    container_image='localhost:30000/refbase:latest'
)

bowtie2_index = ShellTask(
    name="bowtie2-index",
    debug=True,
    script=
    """
    bowtie2-build {inputs.ref} {outputs.idx}/GRCh38_short
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/ref/idx')],
    container_image=bowtie_image_spec
)

bowtie2_align_paired_reads = ShellTask(
    name="bowtie2-align-reads",
    debug=True,
    script=
    """
    bowtie2 -x {inputs.idx} -1 {inputs.read1} -2 {inputs.read2} -S {outputs.sam}
    """,
    inputs=kwtypes(idx=FlyteFile, read1=FlyteFile, read2=FlyteFile),
    output_locs=[OutputLocation(var="sam", var_type=FlyteFile, location='out.sam')],
    container_image=bowtie_image_spec
)

@workflow
def bowtie_wf() -> FlyteDirectory:
    idx = bowtie2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    return idx
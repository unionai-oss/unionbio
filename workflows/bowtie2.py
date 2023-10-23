from flytekit import kwtypes, task, workflow, ImageSpec
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
# from workflows import config, logger

bowtie_image_spec = ImageSpec(
    name="bowtie2",
    apt_packages=["bowtie2"],
    registry="localhost:30000",
    base_image='ghcr.io/pryce-turner/variant-discovery:latest'
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

check_image_spec = ShellTask(
    name="check-image-spec",
    debug=True,
    script=
    """
    which bowtie2
    which fastqc
    """,
    inputs=kwtypes(),
    output_locs=[],
    container_image=bowtie_image_spec
)

# add caching to this task
bowtie2_index = ShellTask(
    name="bowtie2-index",
    debug=True,
    script=
    """
    mkdir {outputs.idx}
    bowtie2-build {inputs.ref} {outputs.idx}/GRCh38_short
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
    container_image=bowtie_image_spec
)

bowtie2_align_paired_reads = ShellTask(
    name="bowtie2-align-reads",
    debug=True,
    script=
    """
    bowtie2 -x {inputs.idx}/GRCh38_short -1 {inputs.read1} -2 {inputs.read2} -S {outputs.sam}
    """,
    inputs=kwtypes(idx=FlyteDirectory, read1=FlyteFile, read2=FlyteFile),
    output_locs=[OutputLocation(var="sam", var_type=FlyteFile, location='out.sam')],
    container_image=bowtie_image_spec
)

@workflow
def bowtie_wf() -> FlyteFile:
    idx = bowtie2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    sam = bowtie2_align_paired_reads(idx=idx, read1='s3://my-s3-bucket/my-data/ERR250683_1.fastq.gz', read2='s3://my-s3-bucket/my-data/ERR250683_2.fastq.gz')
    return sam
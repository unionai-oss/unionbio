from flytekit import kwtypes, task, workflow, ImageSpec
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
# from workflows import config, logger

hisat_image_spec = ImageSpec(
    base_image="refbase:latest",
    apt_packages=["hisat2"]
)

hisat2_index = ShellTask(
    name="hisat2-index",
    debug=True,
    script=
    """
    hisat2-build {inputs.ref} {outputs.idx}/GRCh38_short
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/ref/idx')],
    container_image=hisat_image_spec
)

hisat2_align_paired_reads = ShellTask(
    name="hisat2-align-reads",
    debug=True,
    script=
    """
    hisat2 -x {inputs.idx} -1 {inputs.read1} -2 {inputs.read2} -S {outputs.sam}
    """,
    inputs=kwtypes(idx=FlyteFile, read1=FlyteFile, read2=FlyteFile),
    output_locs=[OutputLocation(var="sam", var_type=FlyteFile, location='out.sam')],
    container_image=hisat_image_spec
)

@workflow
def hisat_wf() -> FlyteDirectory:
    idx = hisat2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    return idx
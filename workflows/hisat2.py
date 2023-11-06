# from flytekit import kwtypes, task, workflow, ImageSpec
# from flytekit.extras.tasks.shell import OutputLocation, ShellTask
# from flytekit.types.file import FlyteFile
# from flytekit.types.directory import FlyteDirectory
# # from workflows import config, logger

# hisat_image_spec = ImageSpec(
#     name="hisat2",
#     apt_packages=["hisat2"],
#     registry="localhost:30000",
#     base_image='ghcr.io/pryce-turner/variant-discovery:latest'
# )

# hisat2_index = ShellTask(
#     name="hisat2-index",
#     debug=True,
#     script=
#     """
#     mkdir {outputs.idx}
#     hisat2-build {inputs.ref} {outputs.idx}/GRCh38_short
#     """,
#     inputs=kwtypes(ref=FlyteFile),
#     output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
#     container_image=hisat_image_spec
# )

# hisat2_align_paired_reads = ShellTask(
#     name="hisat2-align-reads",
#     debug=True,
#     script=
#     """
#     hisat2 -x {inputs.idx}/GRCh38_short -1 {inputs.read1} -2 {inputs.read2} -S {outputs.sam}
#     """,
#     inputs=kwtypes(idx=FlyteDirectory, read1=FlyteFile, read2=FlyteFile),
#     output_locs=[OutputLocation(var="sam", var_type=FlyteFile, location='out.sam')],
#     container_image=hisat_image_spec
# )

# @workflow
# def hisat_wf() -> FlyteFile:
#     idx = hisat2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
#     sam = hisat2_align_paired_reads(idx=idx, read1='s3://my-s3-bucket/my-data/ERR250683_1.fastq.gz', read2='s3://my-s3-bucket/my-data/ERR250683_2.fastq.gz')
#     return sam
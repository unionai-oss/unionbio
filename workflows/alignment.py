from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile

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

bowtie_image_spec = ImageSpec(
    name="bowtie2",
    apt_packages=["bowtie2"],
    registry="localhost:30000",
    base_image='ghcr.io/pryce-turner/variant-discovery:latest'
)

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

hisat_image_spec = ImageSpec(
    name="hisat2",
    apt_packages=["hisat2"],
    registry="localhost:30000",
    base_image='ghcr.io/pryce-turner/variant-discovery:latest'
)

hisat2_index = ShellTask(
    name="hisat2-index",
    debug=True,
    script=
    """
    mkdir {outputs.idx}
    hisat2-build {inputs.ref} {outputs.idx}/GRCh38_short
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
    container_image=hisat_image_spec
)

hisat2_align_paired_reads = ShellTask(
    name="hisat2-align-reads",
    debug=True,
    script=
    """
    hisat2 -x {inputs.idx}/GRCh38_short -1 {inputs.read1} -2 {inputs.read2} -S {outputs.sam}
    """,
    inputs=kwtypes(idx=FlyteDirectory, read1=FlyteFile, read2=FlyteFile),
    output_locs=[OutputLocation(var="sam", var_type=FlyteFile, location='out.sam')],
    container_image=hisat_image_spec
)

multiqc_image_spec = ImageSpec(
    name="multiqc",
    packages=["multiqc"],
    registry="localhost:30000",
    base_image='ghcr.io/pryce-turner/variant-discovery:latest'
)

multiqc = ShellTask(
    name="multiqc",
    debug=True,
    script=
    """
    multiqc {inputs.report_dir} -n {outputs.o}
    """,
    inputs=kwtypes(report_dir=FlyteDirectory),
    output_locs=[OutputLocation(var="o", var_type=FlyteFile, location='/root/multiqc_report.html')],
    container_image=multiqc_image_spec
)

@task(disable_deck=False)
def render_multiqc(report: FlyteFile):
    report_html = open(report, 'r').read()
    current_context().default_deck.append(report_html)

@workflow
def alignment_wf() -> FlyteFile:
    qc = fastqc(seq_dir='s3://my-s3-bucket/my-data/sequences')
    fastp_out = fastp(i1='s3://my-s3-bucket/my-data/sequences/ERR250683_1.fastq.gz', i2='s3://my-s3-bucket/my-data/sequences/ERR250683_2.fastq.gz')
    bowtie2_idx = bowtie2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    bowtie2_sam = bowtie2_align_paired_reads(idx=bowtie2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
    hisat2_idx = hisat2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    hisat2_sam = hisat2_align_paired_reads(idx=hisat2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
    return hisat2_sam
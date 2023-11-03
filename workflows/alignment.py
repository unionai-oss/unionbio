import logging
from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.experimental import eager
from typing import List
from dataclasses import dataclass, asdict
from dataclasses_json import dataclass_json
from pathlib import Path
from flytekit.configuration import Config
from flytekit.remote import FlyteRemote
from mashumaro.mixins.json import DataClassJSONMixin

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s"))
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

base_image = 'localhost:30000/variant-discovery:latest'

@dataclass
class RawSample(DataClassJSONMixin):
    sample: str
    raw_read1: FlyteFile
    raw_read2: FlyteFile

@dataclass
class FiltSample(DataClassJSONMixin):
    sample: str
    filt_read1: FlyteFile
    filt_read2: FlyteFile
    rep: FlyteFile

# fastqc = ShellTask(
#     name="fastqc",
#     debug=True,
#     cache=True,
#     cache_version="1",
#     script=
#     """
#     mkdir {outputs.qc}
#     fastqc {inputs.seq_dir}/*.fastq.gz --outdir={outputs.qc}
#     """,
#     inputs=kwtypes(seq_dir=FlyteDirectory),
#     output_locs=[OutputLocation(var="qc", var_type=FlyteDirectory, location='/root/qc')],
#     container_image=base_image
# )

@task#(cache=True, cache_version="2")
def prepare_samples(seq_dir: FlyteDirectory) -> List[RawSample]:
    samples = {}

    # Fetch FlyteDirectory from object storage and make
    # list of relevant paths
    seq_dir.download()
    all_paths = list(Path(seq_dir.path).rglob('*fastq.gz*'))

    for fp in all_paths:
        
        # Parse paths following 'sample_read.fastq.gz' format
        fn = fp.name
        fullname = fn.split('.')[0]
        sample, mate = fullname.split('_')[0:2]
        
        if not samples.get(sample):
            samples[sample] = RawSample(
                sample=sample,
                raw_read1=FlyteFile(path='/dev/null'),
                raw_read2=FlyteFile(path='/dev/null'),
            )

        print(f'Working on {fn} with mate {mate} for sample {sample}')
        if mate == '1':
            setattr(samples[sample], 'raw_read1', FlyteFile(path=str(fp)))
        elif mate == '2':
            setattr(samples[sample], 'raw_read2', FlyteFile(path=str(fp)))

    return list(samples.values())

fastp = ShellTask(
    name="fastp-shell",
    debug=True,
    # cache=True,
    # cache_version="1",
    container_image=base_image,
    requests=Resources(cpu="1", mem="2Gi"),
    script=
    """
    fastp -i {inputs.i1} -I {inputs.i2} -o {outputs.o1} -O /root/out2.fq.gz
    """,
    inputs=kwtypes(sample=str, i1=FlyteFile, i2=FlyteFile),
    output_locs=[
        OutputLocation(var="o1", var_type=FlyteFile, location='/root/out1.fq.gz'),
        OutputLocation(var="o2", var_type=FlyteFile, location='/root/out2.fq.gz'),
        OutputLocation(var="rep", var_type=FlyteFile, location='/root/fastp.json'),
        ]
)

# bowtie_image_spec = ImageSpec(
#     name="bowtie2",
#     apt_packages=["bowtie2"],
#     registry="localhost:30000",
#     base_image=base_image
# )

# bowtie2_index = ShellTask(
#     name="bowtie2-index",
#     debug=True,
#     script=
#     """
#     mkdir {outputs.idx}
#     bowtie2-build {inputs.ref} {outputs.idx}/GRCh38_short
#     """,
#     inputs=kwtypes(ref=FlyteFile),
#     output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
#     container_image=bowtie_image_spec
# )

# bowtie2_align_paired_reads = ShellTask(
#     name="bowtie2-align-reads",
#     debug=True,
#     script=
#     """
#     bowtie2 -x {inputs.idx}/GRCh38_short -1 {inputs.read1} -2 {inputs.read2} -S {outputs.sam}
#     """,
#     inputs=kwtypes(idx=FlyteDirectory, read1=FlyteFile, read2=FlyteFile),
#     output_locs=[OutputLocation(var="sam", var_type=FlyteFile, location='out.sam')],
#     container_image=bowtie_image_spec
# )

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

# multiqc_image_spec = ImageSpec(
#     name="multiqc",
#     packages=["multiqc"],
#     registry="localhost:30000",
#     base_image='ghcr.io/pryce-turner/variant-discovery:latest'
# )

# multiqc = ShellTask(
#     name="multiqc",
#     debug=True,
#     script=
#     """
#     multiqc {inputs.report_dir} -n {outputs.o}
#     """,
#     inputs=kwtypes(report_dir=FlyteDirectory),
#     output_locs=[OutputLocation(var="o", var_type=FlyteFile, location='/root/multiqc_report.html')],
#     container_image=multiqc_image_spec
# )

# @task(disable_deck=False)
# def render_multiqc(report: FlyteFile):
#     report_html = open(report, 'r').read()
#     current_context().default_deck.append(report_html)

@eager(
    container_image=base_image,
    remote=FlyteRemote(
        config=Config.for_sandbox(),
        default_project="flytesnacks",
        default_domain="development",
    )
)
async def alignment(seq_dir: FlyteDirectory='s3://my-s3-bucket/my-data/single') -> List[FiltSample]:
    # qc = fastqc(seq_dir=seq_dir)
    samples = await prepare_samples(seq_dir=seq_dir)
    filtered_samples = []
    for sample in samples:
        out = await fastp(i1=sample.raw_read1, i2=sample.raw_read2)
    
        filtered_sample = FiltSample(
                sample=sample.sample,
                filt_read1=out.o1,
                filt_read2=out.o2,
                rep=out.rep
            )
        logger.info(f'Created filtered sample with {asdict(filtered_sample)}')
        filtered_samples.append(filtered_sample)
    return filtered_samples

# @dynamic(container_image=base_image)
# def process_samples_bowtie2(samples: List[Sample]):
#     for sample in samples:
        # bowtie2_sam = bowtie2_align_paired_reads(idx=bowtie2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
        # return bowtie2_sam

# @dynamic(container_image=base_image)
# def process_samples_hisat2(samples: List[Sample]):
#     for sample in samples:
        # bowtie2_sam = bowtie2_align_paired_reads(idx=bowtie2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
        # return bowtie2_sam

# @workflow
# def alignment_wf(seq_dir: FlyteDirectory):# -> FlyteFile:
    
    # fastp_out = fastp(i1='s3://my-s3-bucket/my-data/sequences/ERR250683_1.fastq.gz', i2='s3://my-s3-bucket/my-data/sequences/ERR250683_2.fastq.gz')
    # bowtie2_idx = bowtie2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    # bowtie2_alignments = bowtie2_align_paired_reads(idx=bowtie2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
    # hisat2_idx = hisat2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    # hisat2_alignments = hisat2_align_paired_reads(idx=hisat2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
    # return hisat2_sam
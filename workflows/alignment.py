import subprocess
import logging
from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from typing import List
from dataclasses import dataclass, asdict
from dataclasses_json import dataclass_json
from pathlib import Path
from mashumaro.mixins.json import DataClassJSONMixin

# Setup the logger
logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("[%(asctime)s %(levelname)s %(name)s] %(message)s"))
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

base_image = 'ghcr.io/pryce-turner/variant-discovery:latest'

@dataclass
class RawSample(DataClassJSONMixin):
    sample: str
    raw_r1: FlyteFile
    raw_r2: FlyteFile

@dataclass
class FiltSample(DataClassJSONMixin):
    sample: str
    filt_r1: FlyteFile
    filt_r2: FlyteFile
    rep: FlyteFile

@task(container_image=base_image)
def make_filt_sample(indir: FlyteDirectory='s3://my-s3-bucket/my-data/filt-sample') -> FiltSample:
    indir.download()
    print(type(indir.path))
    print(indir.path)
    return FiltSample(
        sample='ERR250683',
        filt_r1=FlyteFile(path=f'{indir.path}/ERR250683_1_filt.fq.gz'),
        filt_r2=FlyteFile(path=f'{indir.path}/ERR250683_2_filt.fq.gz'),
        rep=FlyteFile(path=f'{indir.path}/ERR250683_report.json')
    )

def subproc_raise(command: List[str]) -> str:
    """Execute a command and capture stdout and stderr."""
    try:
        # Execute the command and capture stdout and stderr
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

        # Access the stdout and stderr output
        # stdout_output = result.stdout
        # stderr_output = result.stderr

    except FileNotFoundError as exc:
        print(f"Process failed because the executable could not be found.\n{exc}")
    except subprocess.CalledProcessError as exc:
        print(
            f"Process failed because did not return a successful return code. "
            f"Returned {exc.returncode}\n{exc}"
        )
    except subprocess.TimeoutExpired as exc:
        print(f"Process timed out.\n{exc}")

    
fastqc = ShellTask(
    name="fastqc",
    debug=True,
    cache=True,
    cache_version="1",
    script=
    """
    mkdir {outputs.qc}
    fastqc {inputs.seq_dir}/*.fastq.gz --outdir={outputs.qc}
    """,
    inputs=kwtypes(seq_dir=FlyteDirectory),
    output_locs=[OutputLocation(var="qc", var_type=FlyteDirectory, location='/root/qc')],
    container_image=base_image
)

@task(container_image=base_image)
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
                raw_r1=FlyteFile(path='/dev/null'),
                raw_r2=FlyteFile(path='/dev/null'),
            )

        print(f'Working on {fn} with mate {mate} for sample {sample}')
        if mate == '1':
            setattr(samples[sample], 'raw_r1', FlyteFile(path=str(fp)))
        elif mate == '2':
            setattr(samples[sample], 'raw_r2', FlyteFile(path=str(fp)))

    return list(samples.values())

@task(
    requests=Resources(cpu="1", mem="2Gi"),
    container_image=base_image
    )
def pyfastp(rs: RawSample) -> FiltSample:

    ldir = Path(current_context().working_directory)
    o1p = ldir.joinpath(f'{rs.sample}_1_filt.fq.gz')
    o2p = ldir.joinpath(f'{rs.sample}_2_filt.fq.gz')
    rep = ldir.joinpath(f'{rs.sample}_report.json')

    cmd = [
    "fastp",
    "-i", rs.raw_r1,
    "-I", rs.raw_r2,
    "-o", o1p,
    "-O", o2p,
    "-j", rep
    ]
    
    subproc_raise(cmd)

    return FiltSample(
        sample=rs.sample,
        filt_r1=FlyteFile(path=str(o1p)),
        filt_r2=FlyteFile(path=str(o2p)),
        rep=FlyteFile(path=str(rep))
    )

# this should be map task
@dynamic(container_image=base_image)
def run_fastp(samples: List[RawSample]) -> List[FiltSample]:
    filtered_samples = []
    for sample in samples:
        fs = pyfastp(rs=sample)
        filtered_samples.append(fs)
        logger.info(f'Created filtered sample with {fs}')
    return filtered_samples


bowtie2_image_spec = ImageSpec(
    name="bowtie2",
    apt_packages=["bowtie2"],
    registry="localhost:30000",
    base_image=base_image
)

bowtie2_index = ShellTask(
    name="bowtie2-index",
    debug=True,
    container_image=bowtie2_image_spec,
    script=
    """
    mkdir {outputs.idx}
    bowtie2-build {inputs.ref} {outputs.idx}/GRCh38_short
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
)

@task(container_image=bowtie2_image_spec)
def bowtie2_align_paired_reads(idx: FlyteDirectory, fs: FiltSample) -> FlyteFile:

    idx.download()
    ldir = Path(current_context().working_directory)
    sam = ldir.joinpath(f'{fs.sample}_bowtie2.sam')
    
    cmd = [
        "bowtie2",
        "-x", f"{idx.path}/GRCh38_short",
        "-1", fs.filt_r1,
        "-2", fs.filt_r2,
        "-S", sam
    ]

    subproc_raise(cmd)
    
    return FlyteFile(path=str(sam))

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

# @dynamic(container_image=base_image)
# def compare_aligners(samples: List[Sample]):
#     # map em out
#     for sample in samples:
#         bowtie2_sam = bowtie2_align_paired_reads(idx=bowtie2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
#         hisat2_sam = hisat2_align_paired_reads(idx=bowtie2_idx, read1=fastp_out.o1, read2=fastp_out.o2)
#         return bowtie2_sam, hisat2_sam

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

@dynamic
def bowtie2_align(idx: FlyteDirectory, sample: FiltSample):#samples: List[FiltSample]) -> List[FlyteFile]:
    # sams = []
    # for sample in samples:
    sam = bowtie2_align_paired_reads(idx=idx, fs=sample)
    #     sams.append(sam)
    # return sams

@workflow
def alignment_wf(seq_dir: FlyteDirectory='s3://my-s3-bucket/my-data/single'):# -> FlyteFile:
    # qc = fastqc(seq_dir=seq_dir)
    # samples = prepare_samples(seq_dir=seq_dir)
    # filtered_samples = run_fastp(samples=samples)
    fs = make_filt_sample(indir='s3://my-s3-bucket/my-data/filt-sample')
    bowtie2_idx = bowtie2_index(ref='s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta')
    # bowtie2_align_paired_reads(idx=bowtie2_idx, fs=fs)
    bowtie2_align(idx=bowtie2_idx, sample=fs)
    # hisat2_idx = hisat2_index(ref='s3://my-s3-bucket/my-data/GRCh38_short.fasta')
    # bowtie2_sam, hisat2_sam = compare_aligners(samples=filtered_samples)
    # return hisat2_sam
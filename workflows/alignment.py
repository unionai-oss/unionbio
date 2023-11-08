import os
import subprocess
import logging
import gzip
import hashlib
import shutil
from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.experimental import map_task
from typing import List, Tuple
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

base_image = 'localhost:30000/variant-discovery:latest'
ref_loc = 's3://my-s3-bucket/my-data/refs/GRCh38_short.fasta'
ref_hash = str(hash(ref_loc))[:4]

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
    report: FlyteFile

@dataclass
class SamFile(DataClassJSONMixin):
    sample: str
    sam: FlyteFile
    report: FlyteFile

@task(container_image=base_image)
def make_filt_sample(indir: FlyteDirectory='s3://my-s3-bucket/my-data/filt-sample') -> FiltSample:
    indir.download()
    print(type(indir.path))
    print(indir.path)
    return FiltSample(
        sample='ERR250683',
        filt_r1=FlyteFile(path=f'{indir.path}/ERR250683_1_filt.fq.gz'),
        filt_r2=FlyteFile(path=f'{indir.path}/ERR250683_2_filt.fq.gz'),
        report=FlyteFile(path=f'{indir.path}/ERR250683_report.json')
    )

def subproc_raise(command: List[str]) -> Tuple[str, str]:
    """Execute a command and capture stdout and stderr."""
    try:
        # Execute the command and capture stdout and stderr
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

        # Access the stdout and stderr output
        return result.stdout, result.stderr

    except subprocess.CalledProcessError as e:
        raise Exception(f"Command: {e.cmd}\nFailed with return code {e.returncode}:\n{e.stderr}")
    
    except FileNotFoundError as e:
        raise Exception(f"Process failed because the executable could not be found. \
        Did you specify a container image in the task definition if using custom dependencies?\n{e}")


# human in the loop after fastqc?
# surafce report in flytedeck
# or just have conditional for quality counts before proceeding
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
        report=FlyteFile(path=str(rep))
    )

# this should be map task
@dynamic
def run_fastp(samples: List[RawSample]) -> List[FiltSample]:
    filtered_samples = []
    for sample in samples:
        fs = pyfastp(rs=sample)
        filtered_samples.append(fs)
        logger.info(f'Created filtered sample with {fs}')
    return filtered_samples

bowtie2_index = ShellTask(
    name="bowtie2-index",
    debug=True,
    cache=True,
    cache_version=ref_hash,
    requests=Resources(cpu="4", mem="10Gi"),
    container_image=base_image,
    script=
    """
    mkdir {outputs.idx}
    bowtie2-build {inputs.ref} {outputs.idx}/bt2_idx
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
)

@task(container_image=base_image, requests=Resources(cpu="4", mem="10Gi"))
def bowtie2_align_paired_reads(idx: FlyteDirectory, fs: FiltSample) -> SamFile:

    idx.download()
    ldir = Path(current_context().working_directory)
    sam = ldir.joinpath(f'{fs.sample}_bowtie2.sam')
    rep = ldir.joinpath(f'{fs.sample}_bowtie2_report.txt')
    
    cmd = [
        "bowtie2",
        "-x", f"{idx.path}/bt2_idx",
        "-1", fs.filt_r1,
        "-2", fs.filt_r2,
        "-S", sam
    ]

    stdout, stderr = subproc_raise(cmd)
    
    with open(rep, 'w') as f:
        f.write(stderr)

    return SamFile(
        sample=fs.sample,
        sam=FlyteFile(path=str(sam)),
        report=FlyteFile(path=str(rep))
    )

hisat2_index = ShellTask(
    name="hisat2-index",
    debug=True,
    cache=True,
    cache_version=ref_hash,
    requests=Resources(cpu="4", mem="10Gi"),
    container_image=base_image,
    script=
    """
    mkdir {outputs.idx}
    hisat2-build {inputs.ref} {outputs.idx}/hs2_idx
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteDirectory, location='/root/idx')],
)

@task(container_image=base_image, requests=Resources(cpu="4", mem="10Gi"))
def hisat2_align_paired_reads(idx: FlyteDirectory, fs: FiltSample) -> SamFile:

    idx.download()
    ldir = Path(current_context().working_directory)
    sam = ldir.joinpath(f'{fs.sample}_hisat2.sam')
    rep = ldir.joinpath(f'{fs.sample}_hisat2_report.txt')
    unc_r1 = ldir.joinpath(f'{fs.sample}_1.fq')
    unc_r2 = ldir.joinpath(f'{fs.sample}_2.fq')
    
    with gzip.open(fs.filt_r1, 'rb') as gz_file, open(unc_r1, 'wb') as out_file:
        shutil.copyfileobj(gz_file, out_file)

    with gzip.open(fs.filt_r2, 'rb') as gz_file, open(unc_r2, 'wb') as out_file:
        shutil.copyfileobj(gz_file, out_file)

    cmd = [
        "hisat2",
        "-x", f"{idx.path}/hs2_idx",
        "-1", unc_r1,
        "-2", unc_r2,
        "-S", sam,
        "--summary-file", rep
    ]

    stdout, stderr = subproc_raise(cmd)

    return SamFile(
        sample=fs.sample,
        sam=FlyteFile(path=str(sam)),
        report=FlyteFile(path=str(rep))
    )

@dynamic
def compare_aligners(bt2_idx: FlyteDirectory, hs2_idx: FlyteDirectory, samples: List[FiltSample]) -> List[List[SamFile]]:
    sams = []
    for sample in samples:
        bt2_sam = bowtie2_align_paired_reads(idx=bt2_idx, fs=sample)
        hs2_sam = hisat2_align_paired_reads(idx=hs2_idx, fs=sample)
        pair = [bt2_sam, hs2_sam]
        sams.append(pair)
    return sams

multiqc_image_spec = ImageSpec(
    name="multiqc",
    packages=["multiqc"],
    registry="localhost:30000",
    base_image='ghcr.io/pryce-turner/variant-discovery:latest'
)

@task
def prep_multiqc_ins(fqc: FlyteDirectory, filt_reps: List[FiltSample], sams: List[List[SamFile]]) -> FlyteDirectory:
    # download all the things
    ldir = Path(current_context().working_directory)
    
    fqc.download()
    for f in os.listdir(fqc.path):
        src = os.path.join(fqc.path, f)
        dest = os.path.join(ldir, f)
        shutil.move(src, dest)

    for filt_rep in filt_reps:
        filt_rep.report.download()
        shutil.move(filt_rep.report.path, ldir)

    for pair in sams:
        pair[0].report.download()
        shutil.move(pair[0].report.path, ldir)
        pair[1].report.download()
        shutil.move(pair[1].report.path, ldir)

    return FlyteDirectory(path=str(ldir))

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
def alignment_wf(seq_dir: FlyteDirectory='s3://my-s3-bucket/my-data/sequences'):
    qc = fastqc(seq_dir=seq_dir)
    samples = prepare_samples(seq_dir=seq_dir)
    filtered_samples = map_task(pyfastp)(rs=samples)
    # bowtie2_idx = bowtie2_index(ref=ref_loc)
    # hisat2_idx = hisat2_index(ref=ref_loc)
    # sams = compare_aligners(bt2_idx=bowtie2_idx, hs2_idx=hisat2_idx, samples=filtered_samples)
    # mqc_prep = prep_multiqc_ins(fqc=qc, filt_reps=filtered_samples, sams=sams)
    # mqc = multiqc(report_dir=mqc_prep)
    # render_multiqc(report=mqc)
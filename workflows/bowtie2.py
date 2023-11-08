from pathlib import Path
from flytekit import kwtypes, task, Resources, current_context
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from .config import ref_hash, base_image
from .sample_types import FiltSample, SamFile
from .utils import subproc_raise

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
from pathlib import Path
from flytekit import kwtypes, task, workflow, ImageSpec, Resources, current_context
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from .config import base_image
from .sample_types import FiltSample, RawSample
from .utils import subproc_raise

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
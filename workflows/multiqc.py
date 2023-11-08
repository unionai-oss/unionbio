import os
import shutil
from flytekit import kwtypes, ImageSpec, current_context, task
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from typing import List
from pathlib import Path

from .sample_types import FiltSample, SamFile

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
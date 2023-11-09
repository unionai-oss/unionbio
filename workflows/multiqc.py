import os
import shutil
from flytekit import ImageSpec, current_context, task
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from typing import List
from pathlib import Path

from .sample_types import FiltSample, SamFile
from .utils import subproc_raise

multiqc_image_spec = ImageSpec(
    name="multiqc",
    packages=["multiqc"],
    registry="localhost:30000",
    base_image="ghcr.io/pryce-turner/variant-discovery:latest",
)


@task(container_image=multiqc_image_spec, disable_deck=False)
def render_multiqc(
    fqc: FlyteDirectory, filt_reps: List[FiltSample], sams: List[List[SamFile]]
) -> FlyteFile:
    # download all the things
    ldir = Path(current_context().working_directory)

    fqc.download()
    for f in os.listdir(fqc.path):
        src = ldir.joinpath(fqc.path, f)
        dest = ldir.joinpath(ldir, f)
        shutil.move(src, dest)

    for filt_rep in filt_reps:
        filt_rep.report.download()
        shutil.move(filt_rep.report.path, ldir)

    for pair in sams:
        pair[0].report.download()
        shutil.move(pair[0].report.path, ldir)
        pair[1].report.download()
        shutil.move(pair[1].report.path, ldir)

    final_report = ldir.joinpath("multiqc_report.html")
    mqc_cmd = ["multiqc", str(ldir), "-n", str(final_report)]
    subproc_raise(mqc_cmd)

    report_html = open(final_report, "r").read()
    current_context().default_deck.append(report_html)

    return FlyteFile(path=str(final_report))

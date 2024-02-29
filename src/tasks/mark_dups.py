from typing import List

from flytekit import TaskMetadata, dynamic, kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile

from tasks.sample_types import Alignment
from config import base_image

"""
Identify and remove duplicates from an alignment file using GATK's MarkDuplicates tool.

This function takes in an alignment file, removes the duplicates and writes out
a deduped alignment file.

Args:
    oafn (str): The name of the output deduped alignment file.
    omfn (str): The name of the output deduping metrics file.
    al (FlyteFile): An alignment file containing duplicate reads.

Returns:
    dal (FlyteFile): A deduped alignment file.
    m (FlyteFile): A deduping metrics file.
"""
mark_dups = ShellTask(
    name="mark_dups",
    debug=True,
    metadata=TaskMetadata(retries=3, cache=True, cache_version="1"),
    script="""
    mkdir /tmp/dedup
    "java" \
    "-jar" \
    "/usr/local/bin/gatk" \
    "MarkDuplicates" \
    -I {inputs.al} \
    -O {outputs.dal} \
    -M {outputs.m} \
    """,
    inputs=kwtypes(oafn=str, omfn=str, al=FlyteFile),
    output_locs=[
        OutputLocation(
            var="dal", var_type=FlyteFile, location="/tmp/dedup/{inputs.oafn}"
        ),
        OutputLocation(
            var="m", var_type=FlyteFile, location="/tmp/dedup/{inputs.omfn}"
        ),
    ],
    container_image=base_image,
)


@dynamic
def mark_dups_samples(sams: List[Alignment]) -> List[Alignment]:
    deduped = []
    for i in sams:
        i.deduped = True
        deduped.append(
            mark_dups(
                oafn=i.get_alignment_fname(), omfn=i.get_metrics_fname(), al=i.sam
            )
        )

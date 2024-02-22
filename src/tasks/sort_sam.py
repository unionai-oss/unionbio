from typing import List

from flytekit import TaskMetadata, dynamic, kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile

from tasks.sample_types import SamFile
from config import base_image

# """
# Sort SAM file based on coordinate.

# This function takes a FlyteDirectory object containing raw sequencing data,
# gathers QC metrics using FastQC, and returns a FlyteDirectory object that
# can be crawled with MultiQC to generate a report.

# Args:
#     seq_dir (FlyteDirectory): An S3 prefix containing raw sequencing data to be processed.

# Returns:
#     qc (FlyteDirectory): A directory containing fastqc report output.
# """
sort_sam = ShellTask(
    name="sort_sam",
    debug=True,
    metadata=TaskMetadata(retries=3, cache=True, cache_version="1"),
    script="""
    mkdir /tmp/sort_sam
    "java" \
    "-jar" \
    "/usr/local/bin/gatk" \
    "SortSam" \
    -I {inputs.sam} \
    -O {outputs.o} \
    --SORT_ORDER coordinate \
    """,
    inputs=kwtypes(out_fname=str, sam=FlyteFile),
    output_locs=[
        OutputLocation(
            var="o", var_type=FlyteFile, location="/tmp/sort_sam/{inputs.out_fname}"
        )
    ],
    container_image=base_image,
)


@dynamic
def sort_samples(sams: List[SamFile]) -> List[SamFile]:
    sorted = []
    for i in sams:
        sorted.append(sort_sam(sample=i.sample, sam=i.sam))

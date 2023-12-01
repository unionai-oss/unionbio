from typing import List

from flytekit import TaskMetadata, dynamic, kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile

from .sample_types import SamFile

# """
# Perform quality control using FastQC.

# This function takes a FlyteDirectory object containing raw sequencing data,
# gathers QC metrics using FastQC, and returns a FlyteDirectory object that
# can be crawled with MultiQC to generate a report.

# Args:
#     seq_dir (FlyteDirectory): An S3 prefix containing raw sequencing data to be processed.

# Returns:
#     qc (FlyteDirectory): A directory containing fastqc report output.
# """
mark_dups = ShellTask(
    name="mark_dups",
    debug=True,
    metadata=TaskMetadata(retries=3, cache=True, cache_version="1"),
    script="""
    gatk MarkDuplicatesSpark \
            -I {inputs.sam} \
            -O {outputs.o} \
            --remove-sequencing-duplicates
    """,
    inputs=kwtypes(sample=str, sam=FlyteFile),
    output_locs=[
        OutputLocation(
            var="o", var_type=FlyteFile, location="/root/{inputs.sample}_dedup.sam"
        )
    ],
    container_image="broadinstitute/gatk:4.4.0.0",
)


@dynamic
def mark_dups_samples(sams: List[SamFile]) -> List[SamFile]:
    deduped = []
    for i in sams:
        deduped.append(mark_dups(sample=i.sample, sam=i.sam))

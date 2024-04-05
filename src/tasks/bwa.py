import gzip
import shutil
from pathlib import Path
from flytekit import kwtypes, task, Resources, current_context, TaskMetadata
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from config import ref_hash, base_image, logger

"""
Generate BWA index files from a reference genome.

Args:
    ref (FlyteFile): A FlyteFile object representing the input reference file.

Returns:
    FlyteDirectory: A FlyteDirectory object containing the index files.
"""
bwa_index = ShellTask(
    name="bwa-index",
    debug=True,
    metadata=TaskMetadata(retries=3, cache=True, cache_version=ref_hash),
    requests=Resources(cpu="4", mem="10Gi"),
    container_image=base_image,
    script="""
    mkdir {outputs.idx}
    mv {inputs.ref} {outputs.idx}
    cd {outputs.idx}
    samtools faidx {inputs.ref}
    bwa index {inputs.ref}
    """,
    inputs=kwtypes(ref=FlyteFile),
    output_locs=[
        OutputLocation(var="idx", var_type=FlyteDirectory, location="/tmp/bwa_idx")
    ],
)
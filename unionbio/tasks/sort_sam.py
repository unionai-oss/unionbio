from typing import List

from flytekit import TaskMetadata, dynamic, kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile

from unionbio.datatypes.alignment import Alignment
from unionbio.config import main_img

"""
Sort SAM file based on coordinate.

Args:
    out_fname (str): The name of the output sorted SAM file.
    sam (FlyteFile): An alignment file in SAM format.

Returns:
    o (FlyteFile): A sorted alignment file in SAM format.
"""
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
    container_image=main_img,
)


@dynamic
def sort_samples(sams: List[Alignment]) -> List[Alignment]:
    sorted = []
    for i in sams:
        i.sorted = True
        fname, rep = i.make_filenames()
        sorted.append(sort_sam(out_fname=fname, sam=i.sam))

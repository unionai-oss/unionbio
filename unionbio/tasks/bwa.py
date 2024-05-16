import shutil
from pathlib import Path
from flytekit import kwtypes, task, Resources, current_context, TaskMetadata
from flytekit.extras.tasks.shell import OutputLocation, ShellTask, subproc_execute
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from config import ref_hash, base_image
from datatypes.reference import Reference


@task(
    container_image=base_image,
    requests=Resources(cpu="4", mem="10Gi"),
    cache=True,
    cache_version=ref_hash,
)
def bwa_index(ref_obj: Reference) -> Reference:
    """Indexes a reference genome using BWA.

    Args:
        ref_obj (Reference): The reference object containing the reference genome.

    Returns:
        Reference: The updated reference object with associated index and metadata.
    """
    ref_obj.ref_dir.download()
    sam_idx = [
        'samtools',
        'faidx',
        str(ref_obj.get_ref_path())
    ]
    sam_result = subproc_execute(sam_idx, cwd=ref_obj.ref_dir.path)

    bwa_idx = [
        'bwa',
        'index',
        str(ref_obj.get_ref_path())
    ]
    bwa_result = subproc_execute(bwa_idx, cwd=ref_obj.ref_dir.path)

    setattr(ref_obj, 'index_name', ref_obj.ref_name)
    setattr(ref_obj, 'indexed_with', 'bwa')

    return ref_obj

from flytekit import kwtypes, task, workflow
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile

@task
def fastqc() -> str:
    return "fastqc"

bowtie2_align_reads = ShellTask(
    name="bowtie2-align-reads",
    debug=True,
    script=
    """
    bowtie2 index {inputs.al} -o {outputs.idx}
    """,
    inputs=kwtypes(al=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteFile, location=f"{{inputs.al}}.crai")],
    # container_image=config['current_image']
)

hisat2_align_reads = ShellTask(
    name="hisat2-align-reads",
    debug=True,
    script=
    """
    hisat2 index {inputs.al} -o {outputs.idx}
    """,
    inputs=kwtypes(al=FlyteFile),
    output_locs=[OutputLocation(var="idx", var_type=FlyteFile, location=f"{{inputs.al}}.crai")],
    # container_image=config['current_image']
)

@task
def gather_stats() -> str:
    return "gather_stats"
import os
from union import ImageSpec, task, workflow
from flytekitplugins.slurm import SlurmRemoteScript, SlurmTask, SlurmFunction

image = ImageSpec(
    name="slurm-connector-workflow",
    packages=["git+https://github.com/flyteorg/flytekit.git#plugins/flytekitplugins-slurm"],
    env={"FLYTE_SDK_LOGGING_LEVEL": "10"},
    builder="union",
)

@task(
    container_image=image,
    task_config=SlurmFunction(
        ssh_config={
            "host": "44.223.100.92",
            "username": "ubuntu",
        },
        sbatch_conf={"partition": "debug", "job-name": "tiny-slurm", "output": "/home/ubuntu/fn_task.log"},
        script="""#!/bin/bash -i
echo Run function with sbatch...
# Run the user-defined task function
{task.fn}
""",
    ),
)
def plus_one(x: int) -> int:
    return x + 1


@task(container_image=image)
def greet(year: int) -> str:
    return f"Hello {year}!!!"


@workflow
def fn_wf(x: int) -> str:
    x = plus_one(x=x)
    msg = greet(year=x)
    return msg


local_config = {
    "host": "slurmgpu",
    "username": "ubuntu",
}
remote_config = {
            "host": "44.223.100.92",
            "username": "ubuntu",
        }

bwa = SlurmTask(
    name="bwa",
    task_config=SlurmRemoteScript(
        ssh_config=remote_config,
        batch_script_path="/home/ubuntu/pryce/scripts/bwa.sh",
        batch_script_args=[
            "-r",
            "/home/ubuntu/pryce/outputs/GRCh38_chr21.fasta",
            "-1",
            "/home/ubuntu/pryce/outputs/SRR812824-sub_1.fastq",
            "-2",
            "/home/ubuntu/pryce/outputs/SRR812824-sub_2.fastq",
            "-o",
            "/home/ubuntu/pryce/outputs/SRR812824-sub",
        ],
        sbatch_conf={
            "partition": "debug",
            "job-name": "tiny-slurm",
        }
    )
)

haplocall = SlurmTask(
    name="haplotype-caller",
    task_config=SlurmRemoteScript(
        ssh_config=remote_config,
        batch_script_path="/home/ubuntu/pryce/scripts/haplocaller.sh",
        batch_script_args=[
            "-r",
            "/home/ubuntu/pryce/outputs/GRCh38_chr21.fasta",
            "-b",
            "/home/ubuntu/pryce/outputs/SRR812824-sub.sorted.bam",
            "-o",
            "/home/ubuntu/pryce/outputs/VCFs",
        ],
        sbatch_conf={
            "partition": "debug",
            "job-name": "tiny-slurm",
        }
    )
)

echo_job = SlurmTask(
    name="slurm-task",
    container_image=image,
    task_config=SlurmRemoteScript(
        ssh_config=remote_config,
        batch_script_path="/home/ubuntu/pryce/scripts/echo.sh",
        # batch_script_args=["1"],
        sbatch_conf={
            "partition": "debug",
            "job-name": "tiny-slurm",
        }
    )
)


@workflow
def wf():
    al = bwa()
    call = haplocall()
    al >> call

@workflow
def test_wf():
    echo_job()

if __name__ == "__main__":
    from flytekit.clis.sdk_in_container import pyflyte
    from click.testing import CliRunner

    runner = CliRunner()
    path = os.path.realpath(__file__)

    print(f">>> LOCAL EXEC <<<")
    result = runner.invoke(pyflyte.main, ["run", path, "wf"])
    print(result.output)
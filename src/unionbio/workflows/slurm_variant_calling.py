import os

from flytekit import workflow
from flytekitplugins.slurm import SlurmRemoteScript, SlurmTask

ssh_config = {
    "host": "slurmgpu",
    "username": "ubuntu",
}

echo_1_job = SlurmTask(
    name="slurm-task",
    task_config=SlurmRemoteScript(
        ssh_config=ssh_config,
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

echo_2_job = SlurmTask(
    name="slurm-task",
    task_config=SlurmRemoteScript(
        ssh_config=ssh_config,
        batch_script_path="/home/ubuntu/pryce/scripts/inspect.sh",
        # batch_script_args=["1"],
        sbatch_conf={
            "partition": "debug",
            "job-name": "tiny-slurm",
        }
    )
)


@workflow
def wf() -> None:
    # t1 = echo_1_job()
    t2 = echo_2_job()
    # t1 >> t2

if __name__ == "__main__":
    from flytekit.clis.sdk_in_container import pyflyte
    from click.testing import CliRunner

    runner = CliRunner()
    path = os.path.realpath(__file__)

    print(f">>> LOCAL EXEC <<<")
    result = runner.invoke(pyflyte.main, ["run", path, "wf"])
    print(result.output)
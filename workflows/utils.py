import zipfile
import subprocess
from pathlib import Path
from typing import List, Tuple
from flytekit import task
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.configuration import Config
from flytekit.remote import FlyteRemote

from .config import base_image
from .sample_types import FiltSample, RawSample


@task
def check_fastqc_reports(rep_dir: FlyteDirectory) -> str:
    """
    Check FastQC reports for errors.

    This function checks FastQC reports for errors and raises an exception if any are found.

    Args:
        rep_dir (FlyteDirectory): The input directory containing FastQC reports.
    """
    rep_dir.download()
    all_zips = list(Path(rep_dir.path).rglob("*fastqc.zip*"))

    for p in all_zips:
        with zipfile.ZipFile(p, "r") as zip_file:
            with zip_file.open("summary.txt") as summary:
                contents = summary.read().decode("utf-8")
                if b"FAIL" in contents:
                    return "FAIL"
                elif b"WARN" in contents:
                    return "WARN"

    return "PASS"


def get_remote(local=None, config_file=None):
    """
    Get remote configuration settings and return a remote object.

    This function retrieves remote configuration settings, including the local flag and
    a configuration file, and uses them to create and return a remote object.

    Args:
        local (bool, optional): A flag indicating whether to use local settings. If True,
            the function will use local settings; if False, it will use remote settings.
            Defaults to None, which implies the use of default settings.
        config_file (str, optional): The path to a custom configuration file. If provided,
            this file will be used for configuration settings. Defaults to None.

    Returns:
        Remote: A remote object configured with the specified settings.
    """
    return FlyteRemote(
        config=Config.auto(
            config_file=(
                None
                if local
                else config_file
                if config_file is not None
                else str(Path.home() / ".flyte" / "config-sandbox.yaml")
            )
        ),
        default_project="flytesnacks",
        default_domain="development",
    )


@task(container_image=base_image)
def prepare_samples(seq_dir: FlyteDirectory) -> List[RawSample]:
    """
    Prepare and process raw sequencing data to create a list of RawSample objects.

    This function processes raw sequencing data located in the specified input directory
    and prepares it to create a list of RawSample objects.

    Args:
        seq_dir (FlyteDirectory): The input directory containing raw sequencing data.

    Returns:
        List[RawSample]: A list of RawSample objects representing the processed sequencing data.
    """
    samples = {}

    # Fetch FlyteDirectory from object storage and make
    # list of relevant paths
    seq_dir.download()
    all_paths = list(Path(seq_dir.path).rglob("*fastq.gz*"))

    for fp in all_paths:
        # Parse paths following 'sample_read.fastq.gz' format
        fn = fp.name
        fullname = fn.split(".")[0]
        sample, mate = fullname.split("_")[0:2]

        if not samples.get(sample):
            samples[sample] = RawSample(
                sample=sample,
                raw_r1=FlyteFile(path="/dev/null"),
                raw_r2=FlyteFile(path="/dev/null"),
            )

        print(f"Working on {fn} with mate {mate} for sample {sample}")
        if mate == "1":
            setattr(samples[sample], "raw_r1", FlyteFile(path=str(fp)))
        elif mate == "2":
            setattr(samples[sample], "raw_r2", FlyteFile(path=str(fp)))

    return list(samples.values())


@task(container_image=base_image)
def make_filt_sample(indir: FlyteDirectory) -> FiltSample:
    """
    Create a FiltSample object from input directory.

    This function is used to create a FiltSample object by specifying the input directory
    containing filtered sample data.

    Args:
        indir (FlyteDirectory, optional): The input directory containing filtered sample data.
            Defaults to "s3://my-s3-bucket/my-data/filt-sample".
    """
    indir.download()
    print(type(indir.path))
    print(indir.path)
    return FiltSample(
        sample="ERR250683",
        filt_r1=FlyteFile(path=f"{indir.path}/ERR250683_1_filt.fq.gz"),
        filt_r2=FlyteFile(path=f"{indir.path}/ERR250683_2_filt.fq.gz"),
        report=FlyteFile(path=f"{indir.path}/ERR250683_report.json"),
    )


def subproc_raise(command: List[str]) -> Tuple[str, str]:
    """
    Execute a command and capture its stdout and stderr.
    Args:
        command (List[str]): The command to be executed as a list of strings.
    Returns:
        Tuple[str, str]: A tuple containing the stdout and stderr output of the command.
    Raises:
        Exception: If the command execution fails, this exception is raised with
            details about the command, return code, and stderr output.
        Exception: If the executable is not found, this exception is raised with
            guidance on specifying a container image in the task definition when
            using custom dependencies.
    """
    try:
        # Execute the command and capture stdout and stderr
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

        # Access the stdout and stderr output
        return result.stdout, result.stderr

    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Command: {e.cmd}\nFailed with return code {e.returncode}:\n{e.stderr}"
        )

    except FileNotFoundError as e:
        raise Exception(
            f"""Process failed because the executable could not be found. 
            Did you specify a container image in the task definition if using 
            custom dependencies?\n{e}"""
        )

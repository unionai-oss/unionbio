import subprocess
import logging
from pathlib import Path
from typing import List, Tuple
from flytekit import task
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.configuration import Config
from flytekit.remote import FlyteRemote

from .config import base_image
from .sample_types import FiltSample, RawSample

def get_remote(local=None, config_file=None):
    return FlyteRemote(
        config=Config.auto(
            config_file=(
                None if local
                else config_file if config_file is not None
                else str(Path.home() / ".flyte" / "config-sandbox.yaml")
            )
        ),
        default_project="flytesnacks",
        default_domain="development",
    )

@task(container_image=base_image)
def prepare_samples(seq_dir: FlyteDirectory) -> List[RawSample]:
    samples = {}

    # Fetch FlyteDirectory from object storage and make
    # list of relevant paths
    seq_dir.download()
    all_paths = list(Path(seq_dir.path).rglob('*fastq.gz*'))

    for fp in all_paths:
        
        # Parse paths following 'sample_read.fastq.gz' format
        fn = fp.name
        fullname = fn.split('.')[0]
        sample, mate = fullname.split('_')[0:2]
        
        if not samples.get(sample):
            samples[sample] = RawSample(
                sample=sample,
                raw_r1=FlyteFile(path='/dev/null'),
                raw_r2=FlyteFile(path='/dev/null'),
            )

        print(f'Working on {fn} with mate {mate} for sample {sample}')
        if mate == '1':
            setattr(samples[sample], 'raw_r1', FlyteFile(path=str(fp)))
        elif mate == '2':
            setattr(samples[sample], 'raw_r2', FlyteFile(path=str(fp)))

    return list(samples.values())

@task(container_image=base_image)
def make_filt_sample(indir: FlyteDirectory='s3://my-s3-bucket/my-data/filt-sample') -> FiltSample:
    indir.download()
    print(type(indir.path))
    print(indir.path)
    return FiltSample(
        sample='ERR250683',
        filt_r1=FlyteFile(path=f'{indir.path}/ERR250683_1_filt.fq.gz'),
        filt_r2=FlyteFile(path=f'{indir.path}/ERR250683_2_filt.fq.gz'),
        report=FlyteFile(path=f'{indir.path}/ERR250683_report.json')
    )

def subproc_raise(command: List[str]) -> Tuple[str, str]:
    """Execute a command and capture stdout and stderr."""
    try:
        # Execute the command and capture stdout and stderr
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

        # Access the stdout and stderr output
        return result.stdout, result.stderr

    except subprocess.CalledProcessError as e:
        raise Exception(f"Command: {e.cmd}\nFailed with return code {e.returncode}:\n{e.stderr}")
    
    except FileNotFoundError as e:
        raise Exception(f"Process failed because the executable could not be found. \
        Did you specify a container image in the task definition if using custom dependencies?\n{e}")

import os
import sys
import time
from pathlib import Path
from kubernetes.client.models import (
    V1Container,
    V1PodSpec,
    V1Toleration,
    V1EnvVar,
)
from flytekit import task, workflow, current_context, PodTemplate
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment

sys.path = ['/home/flytekit/workspace/unionbio/src', '/root/micromamba/envs/dev/bin', '/root/unionbio/src/unionbio', '/root', '/root/micromamba/envs/dev/lib/python310.zip', '/root/micromamba/envs/dev/lib/python3.10', '/root/micromamba/envs/dev/lib/python3.10/lib-dynload', '/root/micromamba/envs/dev/lib/python3.10/site-packages']

from unionbio.config import colabfold_img_fqn, logger

DB_LOC = "/root/af_dbs/"
CPU = "14"

actor = ActorEnvironment(
    name="colabfold-actor",
    replica_count=1,
    ttl_seconds=899,
    container_image=colabfold_img_fqn,
)


# @actor.task
@task
def dl_dbs(
    db_uri: str = "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold_dbs_20240927.tar.zst",
    output_loc: str = DB_LOC,
    threads: str = CPU,
) -> str:
    os.makedirs(output_loc, exist_ok=True)

    dl_cmd = [
        "gcloud",
        "storage",
        "cp",
    ]

    if "*" in db_uri:
        dl_cmd += [
            "-r",
            db_uri,
            output_loc,
        ]
    elif db_uri.endswith(".tar.zst"):
        dl_cmd += [
            db_uri,
            "-",
            "|",
            "tar",
            "-I",
            f"'zstd -T{threads}'",
            "--strip-components=1",
            "-x",
            "-C",
            output_loc,
        ]
    else:
        dl_cmd += [
            db_uri,
            output_loc,
        ]

    cmd_str = " ".join(dl_cmd)
    logger.info(f"Downloading databases with command: {cmd_str}")
    start = time.time()
    subproc_execute(command=cmd_str, shell=True)
    logger.info(f"Downloaded in {time.time() - start} seconds")
    logger.debug(f"Database files: {os.listdir(output_loc)}")
    return output_loc


@task
def msa_search(
    seq: FlyteFile = "gs://opta-gcp-dogfood-gcp/bio-assets/P84868.fasta",
):
    from colabfold.mmseqs import search  # type: ignore

    indir = Path(current_context().working_directory).joinpath("inputs")
    outdir = Path(current_context().working_directory).joinpath("outputs")
    os.makedirs(indir, exist_ok=True)
    seq.download()

    # Define the source and destination paths
    source = Path(seq.path)
    destination = indir.joinpath(source.name)

    # Move the file to the destination directory
    source.rename(destination)

    logger.debug(f"Running MMSeqs search on {destination}")
    t = time.time()
    sys.argv = [
        sys.argv[0],
        str(indir),
        str(DB_LOC),
        str(outdir),
        "--db-load-mode",
        "0",
        "--split",
        "0"
    ]
    search.main()

    logger.info(f"Created the following outputs in {time.time() - t} seconds:")
    logger.info(f"MSA files: {os.listdir(outdir)}")

    return FlyteDirectory(path=outdir)


@task
def af_predict(msas: FlyteDirectory) -> FlyteDirectory:
    from colabfold import batch

    outdir = Path(current_context().working_directory).joinpath("outputs")
    msas.download()
    logger.info(f"Running AlphaFold on {os.listdir(msas.path)}")

    sys.argv = [
        sys.argv[0],
        f"--num-models=1",
        msas.path,
        outdir,
    ]
    batch.main()

    return FlyteDirectory(path=outdir)

@workflow
def cf_wf():
    dl = dl_dbs()
    msas = msa_search()
    af = af_predict(msas=msas)

import os
import sys
import time
from pathlib import Path
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
    ttl_seconds=600,
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
    elpased = time.time() - start
    logger.info(f"Downloaded in {elpased} seconds ({elpased/3600} hours)")
    logger.debug(f"Database files: {os.listdir(output_loc)}")
    return output_loc


@task
def cf_search(
    seq: FlyteFile = "gs://opta-gcp-dogfood-gcp/bio-assets/rcsb_pdb_2HHB.fasta",
    db_path: str = DB_LOC,
) -> tuple[FlyteFile, FlyteFile]:

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
    cmd = [
        "colabfold_search",
        "--use-env",
        "1",
        "--use-templates",
        "1",
        "--db-load-mode",
        "2",
        "--db2",
        "pdb100_230517",
        "--threads",
        CPU,
        seq.path,
        db_path,
        str(outdir),
    ]
    subproc_execute(cmd)

    logger.info(f"Created the following outputs in {time.time() - t} seconds:")
    logger.info(f"MSA files: {os.listdir(outdir)}")

    return FlyteDirectory(path=outdir)


@task
def af_predict(hitfile: FlyteFile, msa: FlyteFile) -> FlyteDirectory:

    outdir = Path(current_context().working_directory).joinpath("outputs")
    msa.download()
    hitfile.download()
    logger.info(f"Running AlphaFold on {os.listdir(msas.path)}")

    cmd = [
        "--amber",
        "--templates",
        "--use-gpu-relax",
        "--pdb-hit-file",
        hitfile.path,
        "--local-pdb-path",
        str(Path(DB_LOC).joinpath("pdb/divided")),
        "--random-seed",
        "0",
        msa.path,
        outdir,
    ]
    subproc_execute(cmd)

    return FlyteDirectory(path=outdir)

@workflow
def cf_wf():
    dl = dl_dbs()
    msas = cf_search()
    af = af_predict(msas=msa)

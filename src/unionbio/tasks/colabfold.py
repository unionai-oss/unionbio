import os
import sys
import time
from pathlib import Path
from flytekit import task, workflow, Resources, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import colabfold_img_fqn, logger
from unionbio.types.protein import Protein


DB_LOC = "/root/af_dbs/"
CPU = "15"


actor = ActorEnvironment(
    name="colabfold-actor",
    replica_count=1,
    ttl_seconds=600,
    requests=Resources(
        cpu=CPU,
        mem="32Gi",
        gpu="1",
    ),
    container_image=colabfold_img_fqn,
)

@actor.task
def dl_dbs(
    db_uri: str,
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

@actor.task
def cf_search(seq: FlyteFile, db_loc: str = DB_LOC) -> FlyteDirectory:
    outdir = Path(current_context.working_directory).joinpath("search_out")
    search_cmd = [
        "colabfold_search",
        "--use-env 1",
        "--use-templates 1",
        "--db-load-mode 2",
        "--db2 pdb100_230517",
        "--mmseqs /root/micromamba/envs/dev/bin/mmseqs",
        f"--threads {CPU}",
        seq.path,
        db_loc,
        outdir,
    ]
    cmd_str = " ".join(search_cmd)
    logger.info(f"Running search with command: {cmd_str}")
    start = time.time()
    subproc_execute(command=search_cmd)
    logger.info(f"Search completed in {time.time() - start} seconds")
    logger.debug(f"Search output: {os.listdir(outdir)}")
    return outdir
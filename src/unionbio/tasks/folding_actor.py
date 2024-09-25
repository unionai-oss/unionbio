import os
import sys
import time
from pathlib import Path
from flytekit import task, workflow, Resources
from flytekit.types.file import FlyteFile
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import colabfold_img_fqn, logger
from unionbio.types.protein import Protein


DB_LOC = "/root/af_dbs/"
CPU = "15"


actor = ActorEnvironment(
    name="folding-actor",
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
    db_uri: str = "gs://opta-gcp-dogfood-gcp/bio-assets/af_dbs_test.tar.zst",
    # db_uri: str = "gs://opta-gcp-dogfood-gcp/bio-assets/af_dbs_small.tar.zst"
    output_loc: str = DB_LOC,
    threads: str = CPU,
) -> str:
    os.makedirs(output_loc, exist_ok=True)
    dl_cmd = [
        "gcloud",
        "storage",
        "cp",
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
    cmd_str = " ".join(dl_cmd)
    logger.info(f"Downloading databases with command: {cmd_str}")
    start = time.time()
    subproc_execute(command=cmd_str, shell=True)
    logger.info(f"Downloaded in {time.time() - start} seconds")
    logger.debug(f"Database files: {os.listdir(output_loc)}")
    assert "params" in os.listdir(output_loc)
    return output_loc
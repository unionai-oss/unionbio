import os
import sys
import time
from pathlib import Path
from flytekit import task, workflow, current_context, PodTemplate
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment

# sys.path = ['/home/flytekit/workspace/unionbio/src', '/root/micromamba/envs/dev/bin', '/root/unionbio/src/unionbio', '/root', '/root/micromamba/envs/dev/lib/python310.zip', '/root/micromamba/envs/dev/lib/python3.10', '/root/micromamba/envs/dev/lib/python3.10/lib-dynload', '/root/micromamba/envs/dev/lib/python3.10/site-packages']

from unionbio.config import colabfold_img_fqn, logger

DB_LOC = "/mnt/colabfold"
CPU = "60"

actor = ActorEnvironment(
    name="colabfold-actor",
    replica_count=1,
    ttl_seconds=600,
    container_image=colabfold_img_fqn,
)


# @actor.task
# @task
def gcloud_dl(
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
    elapsed = time.time() - start
    logger.info(f"Downloaded in {elapsed} seconds ({elapsed/3600} hours)")
    logger.debug(f"Database files: {os.listdir(output_loc)}")
    return output_loc

# @task
def s3_sync(
    db_uri: str,
    output_loc: str = DB_LOC,
    retries: int = 5,
) -> str:
    os.makedirs(output_loc, exist_ok=True)

    dl_cmd = [
        "aws",
        "s3",
        "sync",
        db_uri,
        output_loc
    ]

    cmd_str = " ".join(dl_cmd)
    logger.info(f"Downloading databases with command: {cmd_str}")
    start = time.time()

    while retries > 0:
        try:
            subproc_execute(command=cmd_str, shell=True)
        except RuntimeError:
            retries -= 1
            if retries == 0: raise
            continue
    
    elapsed = time.time() - start
    logger.info(f"Downloaded in {elapsed} seconds ({elapsed/3600} hours)")
    logger.debug(f"Database files: {os.listdir(output_loc)}")
    return output_loc


@task
def cf_search(
    seq: FlyteFile,
    db_path: str = DB_LOC,
    outdir: str = None,
) -> tuple[FlyteFile, FlyteFile]:

    outdir = outdir or str(Path(current_context().working_directory).joinpath("outputs"))
    seq.download()

    t = time.time()
    cmd = [
        "colabfold_search",
        "-s",
        "1",
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
        outdir,
    ]
    logger.debug(f"Running MMSeqs search on {seq.path} with command:")
    logger.debug(' '.join(cmd))
    # subproc_execute(cmd)
    logger.info(f"Created the following outputs in {time.time() - t} seconds:")
    logger.info(f"MSA files in {Path(outdir).absolute()}: {os.listdir(outdir)}")

    for fn in os.listdir(outdir):
        path = Path(fn).absolute()
        if path.suffix == ".m8":
            hitfile = FlyteFile(path=str(path))
        if path.suffix == ".a3m":
            msa = FlyteFile(path=str(path))

    logger.debug(f"Returning {hitfile} and {msa}")

    return hitfile, msa


@task
def af_predict(
    hitfile: FlyteFile, 
    msa: FlyteFile, 
    mmcif_loc: str = None,
    outdir: str = None,
    ) -> FlyteDirectory:

    outdir = outdir or str(Path(current_context().working_directory).joinpath("outputs"))
    mmcif_loc = mmcif_loc or str(Path(DB_LOC).joinpath("pdb_mmcif", "mmcif_files"))
    msa.download()
    hitfile.download()
    logger.info(f"Running AlphaFold on {msa.path} and {hitfile.path}")

    cmd = [
        "colabfold_batch",
        "--amber",
        "--templates",
        "--use-gpu-relax",
        "--pdb-hit-file",
        hitfile.path,
        "--local-pdb-path",
        mmcif_loc,
        "--random-seed",
        "0",
        msa.path,
        outdir,
    ]
    logger.debug(f"Executing:")
    logger.debug(' '.join(cmd))
    subproc_execute(cmd)

    return FlyteDirectory(path=outdir)

@workflow
def cf_wf():
    # dl = s3_sync(db_uri="s3://pryce-test/colabfold/")
    hitfile, msa = cf_search(seq="/mnt/folding_io/fastas/P01308.fasta", outdir="/tmp/flyte20mwg4m3/user_spacen0xnt8re/outputs")
    af = af_predict(
        hitfile="/tmp/flyteiwbp4nxz/user_spaceyifdbswm/outputs/sp_P01308_INS_HUMAN_Insulin_OS_Homo_sapiens_OX_9606_GN_INS_PE_1_SV_1_pdb100_230517.m8",
        msa="/tmp/flyteiwbp4nxz/user_spaceyifdbswm/outputs/sp_P01308_INS_HUMAN_Insulin_OS_Homo_sapiens_OX_9606_GN_INS_PE_1_SV_1.a3m",
        mmcif_loc="/mnt/pdb_mmcif/mmcif_files"
        )
    # return af

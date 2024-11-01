import os
import sys
import time
from pathlib import Path
from flytekit import task, workflow, current_context, PodTemplate, Resources, dynamic
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import colabfold_img_fqn, logger

DB_LOC = "/root/colabfold_dbs"
CPU = "30"

actor = ActorEnvironment(
    name="colabfold-actor",
    replica_count=1,
    ttl_seconds=600,
    requests=Resources(
        cpu=CPU,
        mem="100Gi",
        gpu="1",
    ),
    container_image=colabfold_img_fqn,
)


# @task
@actor.task
def gcloud_dl(
    db_uri: str = "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/*",
    output_loc: str = DB_LOC,
) -> str:
    os.makedirs(output_loc, exist_ok=True)

    dl_cmd = [
        "gcloud",
        "storage",
        "rsync",
        "-R",
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
    # return output_loc
    return str(elapsed)

@dynamic
def fetch_all_assets():
    for uri in [
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/cf_envdb/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb100/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/uniref30/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb/",
        # "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb_mmcif/*",
    ]:
        gcloud_dl(db_uri=uri)

@task
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


@actor.task
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
    logger.info(f"MSA files in {Path(outdir).resolve()}: {os.listdir(outdir)}")

    for fn in Path(outdir).iterdir():
        path = fn.resolve()
        if path.suffix == ".m8":
            hitfile = FlyteFile(path=str(path))
        if path.suffix == ".a3m":
            msa = FlyteFile(path=str(path))

    logger.debug(f"Returning {hitfile} and {msa}")

    return hitfile, msa


@actor.task
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

    t = time.time()
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
    logger.info(f"Created the following outputs in {time.time() - t} seconds:")
    logger.info(f"Output files in {Path(outdir).resolve()}: {os.listdir(outdir)}")

    return FlyteDirectory(path=outdir)

@actor.task
def zombie() -> str:
    resurrections = 200
    while True:
        resurrections -= 1
        logger.info("I'm a zombie!")
        time.sleep(30)
        if resurrections == 0:
            return "No more brainsss"
        
@actor.task
def reg_dl(db_uri: FlyteDirectory) -> str:
    s = time.time()
    db_uri.download()
    el = time.time() - s
    return str(el)

@workflow
def cf_wf():# -> FlyteDirectory:
    d1 = gcloud_dl(db_uri = "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb100/") # only one that works
    # d2 = gcloud_dl("gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb/")
    # d4 = gcloud_dl("gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/uniref30/")
    # d1 = gcloud_dl("gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/cf_envdb/")
    # d1 >> d2 >> d3 >> d4
    # fetch_all_assets()
    # z = zombie()
    # hitfile, msa = cf_search(
    #     seq="gs://opta-gcp-dogfood-gcp/bio-assets/P01308.fasta",
    # )
    # af = af_predict(
    #     hitfile=hitfile,
    #     msa=msa,
    # )
    # return af

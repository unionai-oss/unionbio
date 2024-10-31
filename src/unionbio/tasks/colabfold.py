import os
import sys
import time
import shlex
import subprocess
from pathlib import Path
from flytekit import task, workflow, current_context, Resources
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
    db_uri: str = "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold",
    output_loc: str = DB_LOC,
    retries: int = 5,
    prog_interval: int = 60,
) -> str:
    os.makedirs(output_loc, exist_ok=True)
    """
    Syncs files from a GCS source to a destination using gcloud storage rsync with retries and progress logger.

    Args:
        source (str): Source GCS path (e.g., 'gs://your-bucket/path').
        destination (str): Local destination path or another GCS path.
        retries (int): Number of retry attempts on failure.
    """

    # Build the gcloud rsync command
    # Calculate the total size of objects to download
    logger.debug("Calculating size of objects to download...")
    gsutil_size_command = f"gsutil du -s {db_uri}"
    # total_size_output = subprocess.check_output(shlex.split(gsutil_size_command))
    # total_size_bytes = int(total_size_output.split()[0])
    # logger.info(f"Total size to download: {total_size_bytes} bytes or {total_size_bytes / (1024 ** 3):.2f} GB")
    
    command = ["gcloud", "storage", "rsync", "-r", db_uri, output_loc]
    
    attempt = 0
    while attempt < retries:
        attempt += 1
        logger.info(f"Starting rsync attempt {attempt} of {retries}")
        
        # Run the rsync command with progress
        logger.debug("Running rsync command:")
        logger.debug(' '.join(command))
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        while process.poll() is None:  # While the download is still running
            # Calculate the current size of downloaded files
            current_size_output = subprocess.check_output(shlex.split(f"du -sb {output_loc}"))
            current_size_bytes = int(current_size_output.split()[0])

            # Calculate and print the download progress
            # progress = (current_size_bytes / total_size_bytes) * 100
            # print(f"Downloaded {current_size_bytes / (1024 ** 3):.2f} GB ({progress:.2f}%)")

            # Wait for the specified interval before checking again
            time.sleep(prog_interval)
            
            # Check if rsync succeeded
            process.poll()
            logger.debug(f"Return code: {process.returncode}")
            if process.returncode is None:
                continue
            elif process.returncode == 0:
                logger.info("Rsync completed successfully.")
                process.communicate()
                return output_loc
            else:
                # Log the error and retry if failed
                logger.error(f"Rsync attempt {attempt} failed with error:\n{process.stderr.read()}")
                if attempt < retries:
                    logger.info(f"Retrying in 5 seconds...")
                    time.sleep(5)
                else:
                    logger.error("Max retries reached. Rsync failed.")
                    process.terminate()

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
    outdir: str | None = None,
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
    mmcif_loc: str | None = None,
    outdir: str | None = None,
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

@workflow
def cf_wf() -> FlyteDirectory:
    dl = gcloud_dl()
    hitfile, msa = cf_search(
        seq="gs://opta-gcp-dogfood-gcp/bio-assets/P01308.fasta",
        db_path=dl,
    )
    af = af_predict(
        hitfile=hitfile,
        msa=msa,
    )
    return af

# @actor.task
# def debug() -> str:
#     from unionbio.config import logger
#     return str(sys.path)
import os
import sys
import time
from pathlib import Path
from typing import Tuple, List
from flytekit import task, workflow, Resources, current_context
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import colabfold_img_fqn, logger
from unionbio.types.protein import Protein


DB_LOC = "/root/af_dbs/"
CPU = "15"
os.environ["ALPHAFOLD_DATA_DIR"] = DB_LOC


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
def generate_msas(
    seq: FlyteFile, db_loc: str = DB_LOC
) -> tuple[FlyteDirectory, list[str], dict[str, str]]:
    from colabfold.batch import get_msa_and_templates

    outdir = Path(current_context.working_directory).joinpath("search_out")
    start = time.time()
    a3m_files, template_results = get_msa_and_templates(
        seq.path,
        msa_mode="mmseqs2_uniref_env",
        use_templates=True,
        msa_output_dir=outdir,
    )
    logger.info(f"Search completed in {time.time() - start} seconds")
    logger.debug(f"Search output: {os.listdir(outdir)}")
    return FlyteDirectory(outdir), a3m_files, template_results


@actor.task
def predict_structure(
    seq: FlyteFile,
    msa: FlyteDirectory,
    a3m_files: list[str],
    templates: dict[str, str],
) -> FlyteDirectory:
    from colabfold.predict import run_alphafold

    # Check if MSA files are present e.g. Actor mode
    if not os.listdir(msa.path):
        msa.download()

    results = Path(current_context.working_directory).joinpath("af_results")
    start = time.time()
    prediction_results = run_alphafold(
        fasta_path=seq.path,
        a3m_files=a3m_files,
        templates_results=templates,
        result_dir=results,
        model_type="alphafold2_multimer_v3",
    )
    logger.info(f"Prediction completed in {time.time() - start} seconds")
    return FlyteDirectory(path=results)

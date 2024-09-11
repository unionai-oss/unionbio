import os
from flytekit import task, workflow, Resources
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import alphafold_img_fqn

actor = ActorEnvironment(
    name="af-actor",
    replica_count=1,
    parallelism=1,
    backlog_length=10,
    ttl_seconds=600,
    requests=Resources(
        cpu="8",
        mem="8Gi",
        ephemeral_storage="600Gi",
        gpu="1",
    ),
    container_image=alphafold_img_fqn,
)


@actor.task
def dl_dbs(
    script_loc: str = "/app/alphafold/scripts/download_all_data.sh",
    output_loc: str = "/mnt/af_dbs",
    reduced: bool = False,
) -> str:
    os.makedirs(output_loc, exist_ok=True)
    dl_cmd = [
        script_loc,
        output_loc,
    ]
    if reduced:
        dl_cmd.append("reduced_dbs")
    subproc_execute(command=dl_cmd)
    return "\n".join(os.listdir(output_loc))


@actor.task
def run_af(entry: str = "/app/run_alphafold.sh"):
    cfg = {
        "use_gpu": True,  # 'Enable NVIDIA runtime to run with GPUs.'
        "models_to_relax": "best",  # ['best', 'all', 'none'],
        # 'The models to run the final relaxation step on. '
        # 'If `all`, all models are relaxed, which may be time '
        # 'consuming. If `best`, only the most confident model is '
        # 'relaxed. If `none`, relaxation is not run. Turning off '
        # 'relaxation might result in predictions with '
        # 'distracting stereochemical violations but might help '
        # 'in case you are having issues with the relaxation '
        # 'stage.'
        "enable_gpu_relax": True,  #'Run relax on GPU if GPU is enabled.'
        "gpu_devices": "all",
        # 'Comma separated list of devices to pass to NVIDIA_VISIBLE_DEVICES.'
        "fasta_paths": None,
        # 'Paths to FASTA files, each containing a prediction '
        # 'target that will be folded one after another. If a FASTA file contains '
        # 'multiple sequences, then it will be folded as a multimer. Paths should be '
        # 'separated by commas. All FASTA paths must have a unique basename as the '
        # 'basename is used to name the output directories for each prediction.'
        "output_dir": "/tmp/alphafold",
        # 'Path to a directory that will store the results.'
        "data_dir": None,
        # 'Path to directory with supporting data: AlphaFold parameters and genetic '
        # 'and template databases. Set to the target of download_all_databases.sh.'
        "max_template_date": None,
        # 'Maximum template release date to consider (ISO-8601 format: YYYY-MM-DD). '
        # 'Important if folding historical test sets.'
        "db_preset": "full_dbs",  # ['full_dbs', 'reduced_dbs'],
        # 'Choose preset MSA database configuration - smaller genetic database '
        # 'config (reduced_dbs) or full genetic database config (full_dbs)'
        "model_preset": "monomer",
        # ['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer'],
        # 'Choose preset model configuration - the monomer model, the monomer model '
        # 'with extra ensembling, monomer model with pTM head, or multimer model'
        "num_multimer_predictions_per_model": 5,
        # 'How many predictions (each with a different random seed) will be '
        # 'generated per model. E.g. if this is 2 and there are 5 '
        # 'models then there will be 10 predictions per input. '
        # 'Note: this FLAG only applies if model_preset=multimer'
        "benchmark": False,
        # 'Run multiple JAX model evaluations to obtain a timing that excludes the '
        # 'compilation time, which should be more indicative of the time required '
        # 'for inferencing many proteins.'
        "use_precomputed_msas": False,
        # 'Whether to read MSAs that have been written to disk instead of running '
        # 'the MSA tools. The MSA files are looked up in the output directory, so it '
        # 'must stay the same between multiple runs that are to reuse the MSAs. '
        # 'WARNING: This will not check if the sequence, database or configuration '
        # 'have changed.'
    }

    af_cmd = [
        entry,
    ]

import os
import time
from pathlib import Path
from flytekit import task, workflow, Resources
from flytekit.types.file import FlyteFile
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import alphafold_img_fqn
from unionbio.types.protein import Protein

actor = ActorEnvironment(
    name="af-actor",
    replica_count=1,
    parallelism=1,
    backlog_length=10,
    ttl_seconds=600,
    requests=Resources(
        cpu="8",
        mem="8Gi",
        ephemeral_storage="300Gi",
        gpu="1",
    ),
    container_image=alphafold_img_fqn,
)

DB_LOC = "/mnt/af_dbs"

@actor.task
# @task
def dl_dbs(
    script_loc: str = "/app/alphafold/scripts/download_all_data.sh",
    output_loc: str = DB_LOC,
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
# @task
def run_af(
    # fastas: list[Protein],
    fasta: FlyteFile="s3://union-cloud-oc-staging-dogfood/bio-assets/P68871_sequence.fasta",
    entry: str = "/app/run_alphafold.sh",
    db_dir: str = DB_LOC,
    cfg_ov: dict | None = None,
    db_cfg_ov: dict | None = None,
):
    def_env = {
        'NVIDIA_VISIBLE_DEVICES': "all", # 'Comma separated list of devices to pass to NVIDIA_VISIBLE_DEVICES.'
        # The following flags allow us to make predictions on proteins that
        # would typically be too long to fit into GPU memory.
        'TF_FORCE_UNIFIED_MEMORY': '1',
        'XLA_PYTHON_CLIENT_MEM_FRACTION': '4.0',
    }
    for var, val in def_env.items():
        if var not in os.environ:
            os.environ[var] = val

    cfg = {
        "use_gpu": True,  # 'Enable NVIDIA runtime to run with GPUs.'
        "enable_gpu_relax": True,  #'Run relax on GPU if GPU is enabled.'
        "models_to_relax": "best",  # ['best', 'all', 'none'],
        # 'The models to run the final relaxation step on. '
        # 'If `all`, all models are relaxed, which may be time '
        # 'consuming. If `best`, only the most confident model is '
        # 'relaxed. If `none`, relaxation is not run. Turning off '
        # 'relaxation might result in predictions with '
        # 'distracting stereochemical violations but might help '
        # 'in case you are having issues with the relaxation '
        # 'stage.'
        "output_dir": "/tmp/alphafold",
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
    cfg = {**cfg, **cfg_ov}

    # Set up the GPU relaxation flag and remove it's components from cfg dict
    cfg["use_gpu_relax"] = cfg["use_gpu"] and cfg["enable_gpu_relax"]
    del cfg["enable_gpu_relax"]
    del cfg["use_gpu"]

    # You can individually override the following paths if you have placed the
    # data in locations other than the default dir
    db_cfg = {
        "data_dir": db_dir,
        # Path to the Uniref90 database for use by JackHMMER.
        "uniref90_database_path": db_dir.joinpath("uniref90", "uniref90.fasta"),
        # Path to the MGnify database for use by JackHMMER.
        "mgnify_database_path": db_dir.joinpath("mgnify", "mgy_clusters_2022_05.fa"),
        # Path to a directory with template mmCIF structures, each named <pdb_id>.cif.
        "template_mmcif_dir": db_dir.joinpath("pdb_mmcif", "mmcif_fi: db_dir.les"),
        # Path to a file mapping obsolete PDB IDs to their replacements.
        "obsolete_pdbs_path": db_dir.joinpath("pdb_mmcif", "obsolete.: db_dir.dat"),
    }

    if db_cfg.get("model_preset") == "multimer":
        # Path to the PDB seqres database for use by hmmsearch.
        db_cfg["pdb_seqres_database_path"] = db_dir.joinpath("pdb_seqres", "pdb_seqres.txt"),
        # Path to the Uniprot database for use by JackHMMER.
        db_cfg["uniprot_database_path"] = db_dir.joinpath("uniprot", "uniprot.fasta"),
    else:
        # Path to the PDB70 database for use by HHsearch.
        db_cfg["pdb70_database_path"] = db_dir.joinpath("pdb70", "pd: db_dir.b70"),

    if cfg.get("db_preset") == "reduced_dbs":
        # Path to the Small BFD database for use by JackHMMER.
        db_cfg["small_bfd_database_path"] = db_dir.joinpath(
            "small_bfd", "bfd-first_non_consensus_sequences.fasta"
        ),
    else:
        # Path to the Uniref30 database for use by HHblits.
        db_cfg["uniref30_database_path"] = db_dir.joinpath("uniref30", "UniRef30_2021_03"),
        # Path to the BFD database for use by HHblits.
        db_cfg["bfd_database_path"] = db_dir.joinpath(
            "bfd", "bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
        ),
    db_cfg = {**db_cfg, **db_cfg_ov}

    af_cmd = [
        entry,
        "--logstostderr",
    ]

    # fasta_str = ",".join([str(fasta.sequence_fasta.path) for fasta in fastas])
    fasta_str = str(fasta.path)
    af_cmd.append(f"--fasta_paths={fasta_str}")
    # 'Paths to FASTA files, each containing a prediction '
    # 'target that will be folded one after another. If a FASTA file contains '
    # 'multiple sequences, then it will be folded as a multimer. Paths should be '
    # 'separated by commas. All FASTA paths must have a unique basename as the '
    # 'basename is used to name the output directories for each prediction.'

    for k, v in db_cfg.items():
        af_cmd.append(f"--{k}={v}")

    for k, v in cfg.items():
        af_cmd.append(f"--{k}={v}")

    cmd_str = " ".join(af_cmd)
    print(cmd_str)

    time.sleep(3600)

@workflow
def af_wf():
    dl_dbs(reduced=True)
    run_af(
        fasta="s3://union-cloud-oc-staging-dogfood/bio-assets/P68871_sequence.fasta",
        cfg_ov={"db_preset": "reduced_dbs"}
        )
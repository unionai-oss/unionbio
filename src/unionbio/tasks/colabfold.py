import os
import time
import base64
from io import BytesIO
from pathlib import Path
from flytekit import task, workflow, current_context, Resources, Deck
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

@actor.task
def sync_dbs(
    uris: list[str] = [
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/cf_envdb/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb100/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/uniref30/",
    ],
    output_loc: str = DB_LOC,
) -> str:
    
    os.makedirs(output_loc, exist_ok=True)
    
    for uri in uris:
        dl_cmd = [
            "gcloud",
            "storage",
            "rsync",
            "-R",
            uri,
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


@actor.task
def sync_mmcif(
    uri: str = "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/mmcif_tar/",
    output_loc: str | None = None,
) -> str:
    
    import zstandard as zstd
    import tarfile
    
    output_loc = output_loc or str(Path(DB_LOC).joinpath("pdb"))
    temp_dir = "/tmp/mmcif/"
    os.makedirs(output_loc, exist_ok=True)

    dl_cmd = [
        "gcloud",
        "storage",
        "rsync",
        "-R",
        uri,
        temp_dir,
    ]
            
    cmd_str = " ".join(dl_cmd)
    logger.info(f"Downloading databases with command: {cmd_str}")
    start = time.time()
    subproc_execute(command=cmd_str, shell=True)
    elapsed = time.time() - start
    logger.info(f"Downloaded in {elapsed} seconds ({elapsed/3600} hours)")
    logger.debug(f"Database files: {os.listdir(output_loc)}")

    # Loop over each file in the source directory
    for filename in os.listdir(temp_dir):
        if filename.endswith(".tar.zst"):
            src_path = os.path.join(temp_dir, filename)
            print(f"Decompressing {src_path}...")

            # Decompress the .zst file
            with open(src_path, 'rb') as compressed_file:
                decompressor = zstd.ZstdDecompressor()
                with decompressor.stream_reader(compressed_file) as decompressed_stream:
                    # Extract the .tar content from the decompressed stream
                    with tarfile.open(fileobj=decompressed_stream, mode='r:') as tar:
                        tar.extractall(path=output_loc)
            
            print(f"Decompressed and extracted {src_path} to {output_loc}")

    return output_loc

@actor.task
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
    subproc_execute(cmd)
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


@task(enable_deck=True, container_image=colabfold_img_fqn)
def visualize(af_res: FlyteDirectory) -> FlyteFile:
    
    import plotly
    from graphein.protein.config import ProteinGraphConfig 
    from graphein.protein.graphs import construct_graph
    from graphein.protein.visualisation import plotly_protein_structure_graph
    from PIL import Image
    
    af_res.download()

    # Select the highest confidence relaxed model
    pdb = list(Path(af_res.path).glob("*relaxed_rank_001*"))[0]
    config = ProteinGraphConfig()
    g = construct_graph(config=config, path=pdb)
    p = plotly_protein_structure_graph(
        g,
        colour_edges_by="kind",
        colour_nodes_by="seq_position",
        label_node_ids=False,
        plot_title="Peptide backbone graph. Nodes coloured by position.",
        node_size_multiplier=1
        )
    
    # Append structure to deck
    prot_deck = Deck("Structure")
    html = plotly.io.to_html(p)
    prot_deck.append(html)

    # Append prediction quality images to deck
    pred_qual = Deck("Prediction Quality")
    for i in list(Path(af_res.path).glob("*.png")):
        with Image.open(i) as img:
            img_bytes = BytesIO()
            img.save(img_bytes, format="PNG")
            image_base64 = base64.b64encode(img_bytes.getvalue()).decode()
            pred_qual.append(f'<html><body><img src="data:image/png;base64,{image_base64}"/></body></html>')

    # Write structure to file
    of = "structure.html"
    with open(of, "w") as f:
        p.write_html(f)

    return FlyteFile(path=of)
        

@workflow
def cf_wf() -> FlyteFile:
    d1 = gcloud_dl(db_uri = "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb100/") # only one that works
    d2 = gcloud_dl("gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb/")
    d4 = gcloud_dl("gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/uniref30/")
    d1 = gcloud_dl("gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/cf_envdb/")
    hitfile, msa = cf_search(
        seq="gs://opta-gcp-dogfood-gcp/bio-assets/P01308.fasta",
    )
    af = af_predict(
        hitfile=hitfile,
        msa=msa,
    )
    plot = visualize(af_res="gs://opta-gcp-dogfood-gcp/bio-assets/hemoglobin/fold_out/")
    return plot

# @actor.task
# def debug() -> str:
#     from unionbio.config import logger
#     return str(sys.path)
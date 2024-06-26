# Ok, we are ready to begin and import all tools we will use in this colab:
import biotite.structure.io as bsio
import torch
from Bio import SeqIO
from flytekit import Resources, task
from flytekit.extras.tasks.shell import subproc_execute
from flytekit.types.file import FlyteFile
from transformers import AutoTokenizer, EsmForProteinFolding
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.protein import Protein
from unionbio.config import logger, folding_img_fqn


@task(container_image=folding_img_fqn)
def test_unionbio_install() -> str:
    return "UnionBio package installed successfully."


@task(container_image=folding_img_fqn)
def prodigal_predict(in_seq: Reads) -> Protein:
    """
    Predicts protein sequences from a DNA sequence using Prodigal.
    """
    seq = in_seq.uread
    seq.download()

    prot = Protein(name=in_seq.sample)
    prot_out = prot.get_prot_fname()
    genes_out = prot.get_genes_fname()
    logger.debug(f"Initialized protein object: {prot}")

    cmd = [
        "prodigal",
        "-i",
        seq.path,
        "-a",
        prot_out,
        "-o",
        genes_out,
        "-p",
        "meta",
        "-f",
        "gff",
    ]
    logger.debug(f"Running Prodigal: {' '.join(cmd)}")
    result = subproc_execute(cmd)
    logger.info(
        f"Command exited with code {result.returncode} and stdout {result.output}"
    )

    setattr(prot, "protein", prot_out)
    setattr(prot, "genes", genes_out)
    logger.info(f"Returning protein object: {prot}")

    return prot


@task(container_image=folding_img_fqn, requests=Resources(gpu="1"))
def esm_fold(prot: FlyteFile) -> FlyteFile:
    esmfold = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1",
        low_cpu_mem_usage=True,  # we set this flag to save some RAM during loading
    )
    # esmfold = esmfold.to(device) # -> move ESMFold to the GPU
    esmfold.esm = esmfold.esm.half()  # -> make sure we use a lightweight precision
    esmfold_tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")

    # 1. Select a predicted protein-coding gene that we want to fold:
    i = 0  # lets just use the first one.
    protein_record = list(SeqIO.parse(prot.path, "fasta"))[i]
    protein_seq = str(protein_record.seq)[:-1]  # remove stop codon

    # 2. Tokenize the gene to make it digestible for ESMFold:
    esmfold_in = esmfold_tokenizer(
        [protein_seq], return_tensors="pt", add_special_tokens=False
    )

    # 3. Feed the tokenized sequence to ESMfold:
    # (in PyTorch's inference_mode to avoid any costly gradient computations)
    with torch.inference_mode():
        esmfold_out = esmfold(**esmfold_in.to("cuda"))
        esmfold_out_pdb = esmfold.output_to_pdb(esmfold_out)[0]

    protein_structure_pdb = "protein_structure.pdb"
    with open(protein_structure_pdb, "w") as f:
        f.write(esmfold_out_pdb)

    protein_structure = bsio.load_structure(
        protein_structure_pdb, extra_fields=["b_factor"]
    )
    print("Generated protein structure: ", protein_structure)

    return FlyteFile(path=protein_structure_pdb)

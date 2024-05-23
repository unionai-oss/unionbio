# Ok, we are ready to begin and import all tools we will use in this colab:
import os
import numpy as np
import matplotlib.pyplot as plt
import py3Dmol
import biotite.structure.io as bsio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from flytekit import ImageSpec, Resources, task, workflow
from flytekit.extras.tasks.shell import subproc_execute
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.protein import Protein
from unionbio.config import folding_img, logger

@task(container_image=folding_img)
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
        "-i", seq.path,
        "-a", prot_out,
        "-o", genes_out,
        "-p", "meta",
        "-f", "gff",
    ]
    logger.debug(f"Running Prodigal: {' '.join(cmd)}")
    result = subproc_execute(cmd)
    logger.info(f"Command exited with code {result.returncode} and stdout {result.output}")

    setattr(prot, "protein", prot_out)
    setattr(prot, "genes", genes_out)
    logger.info(f"Returning protein object: {prot}")

    return prot

# Ok, we are ready to begin and import all tools we will use in this colab:
import os
import numpy as np
import matplotlib.pyplot as plt
import py3Dmol # -> used to visualize protein structures
import biotite.structure.io as bsio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from flytekit import ImageSpec, Resources, task, workflow
from datatypes.reads import Reads
from datatypes.protein import Protein
from config import test_assets

folding_img = ImageSpec(
    image_name="protein_folding",
    base_image="ghcr.io/flyteorg/flytekit:py3.11-1.12.0",
    packages=["biopython", "biotite", "transformers", "py3Dmol"],
    conda_channels=["bioconda"],
    conda_packages=["prodigal"],
    registry="ghcr.io/pryce-turner"
)

@task(container_image=folding_img)
def prodigal_predict(in_seq: Reads) -> Protein:
    """
    Predicts protein sequences from a DNA sequence using Prodigal.
    """
    seq = in_seq.read1
    seq.download()

    prot = Protein(name=in_seq.sample)
    prot_out = prot.get_prot_fname()
    genes_out = prot.get_genes_fname()

    cmd = [
        "prodigal",
        "-i", seq.path,
        "-a", prot_out,
        "-o", genes_out,
    ]

    setattr(prot, "protein", prot_out)
    setattr(prot, "genes", genes_out)

    return prot

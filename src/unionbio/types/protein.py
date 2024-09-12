from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile


@dataclass
class Protein(DataClassJSONMixin):
    """
    Represents a protein sequence and its associated metadata.

    This class defines the structure for representing a protein 
    sequence along with attributes that describe the associated 
    metadata.

    Attributes:
        name (str): The name or identifier of the protein sequence.
        sequence_fasta (FlyteFile): A FlyteFile object representing 
            the path to the amino acid sequence fasta file.
        genes (FlyteFile): A FlyteFile object representing the path 
            to the gene sequence file.
    """

    name: str
    sequence_fasta: FlyteFile | None = None
    genes: FlyteFile | None = None

    def get_prot_fname(self):
        return f"{self.name}_sequence.fasta"

    def get_genes_fname(self):
        return f"{self.name}_genes.gff"

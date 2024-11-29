from pathlib import Path
from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from unionbio.config import logger
from unionbio.tasks.helpers import filter_dir

@dataclass
class Protein(DataClassJSONMixin):
    """
    Represents a protein sequence and its associated metadata.

    This class defines the structure for representing a protein sequence along with attributes
    that describe the associated metadata.

    Attributes:
        name (str): The name or identifier of the protein sequence.
        protein (FlyteFile): A FlyteFile object representing the path to the protein sequence file.
        genes (FlyteFile): A FlyteFile object representing the path to the gene sequence file.
    """

    sample: str
    sequence: FlyteFile | None = None
    msa: FlyteFile | None = None
    hitfile: FlyteFile | None = None
    genes: FlyteFile | None = None

    def get_prot_fname(self):
        return f"{self.name}.fasta"

    def get_genes_fname(self):
        return f"{self.name}.gff"

    @classmethod
    def make_all(
        cls,
        dir: Path,
        include: list[str] = ["*.fasta*", "*.gff", "*.m8", "*.a3m"],
        exclude: list[str] = [],
    ) -> list:
        samples = {}
        for fp in filter_dir(dir, include=include, exclude=exclude):
            sample = fp.stem

            if sample not in samples:
                samples[sample] = cls(sample=sample)

            if fp.suffix == ".fasta":
                samples[sample].sequence = FlyteFile(path=str(fp))
            elif fp.suffix == ".m8":
                samples[sample].alignment_idx = FlyteFile(path=str(fp))

        logger.info(f"Created following Alignment objects from {dir}: {samples}")
        return list(samples.values())

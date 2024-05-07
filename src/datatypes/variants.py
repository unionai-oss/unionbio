from pathlib import Path
from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from config import logger


@dataclass
class VCF(DataClassJSONMixin):
    """
    Represents a VCF (Variant Call Format) file and its associated sample and index.

    This class defines the structure for representing a VCF file along with attributes
    that describe the associated sample and index.

    Attributes:
        sample (str): The name or identifier of the sample to which the VCF file belongs.
        caller (str): The name of the variant caller used to generate the VCF.
        vcf (FlyteFile): The VCF file.
        vcf_idx (FlyteFile): The index file for the VCF.

    """

    sample: str
    caller: str
    vcf: FlyteFile | None = None
    vcf_idx: FlyteFile | None = None

    def _get_state_str(self):
        state = f"{self.sample}_{self.caller}"
        return state

    def get_vcf_fname(self):
        return f"{self._get_state_str()}.vcf"

    def get_vcf_idx_fname(self):
        return f"{self._get_state_str()}.vcf.tbi"

    def dl_all(self):
        self.vcf.download()
        self.vcf_idx.download()

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        pattern = "*vcf*"
        dir_contents = list(dir.rglob(pattern))
        logger.info(f"Found following VCF files in {dir} matching {pattern}: {dir_contents}")
        for fp in dir_contents:
            sample, aligner = fp.stem.split("_")[0:2]

            if sample not in samples:
                samples[sample] = VCF(sample=sample, aligner=aligner)

            if "vcf" in fp.name:
                setattr(samples[sample], "vcf", FlyteFile(path=str(fp)))
            elif "tbi" in fp.name:
                setattr(samples[sample], "vcf_idx", FlyteFile(path=str(fp)))

        logger.info(f"Created following VCF objects from {dir}: {samples}")
        return list(samples.values())

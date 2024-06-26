from pathlib import Path
from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from unionbio.config import logger


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
        return f"{self._get_state_str()}.vcf.gz"

    def get_vcf_idx_fname(self):
        return f"{self._get_state_str()}.vcf.gz.tbi"

    def dl_all(self, workdir: Path):
        v_loc = Path(self.vcf.download())
        i_loc = Path(self.vcf_idx.download())
        v_name = v_loc.name
        i_name = i_loc.name
        v_new = workdir.joinpath(v_name)
        i_new = workdir.joinpath(i_name)
        v_loc.rename(v_new)
        i_loc.rename(i_new)
        self.vcf = FlyteFile(path=str(v_new))
        self.vcf_idx = FlyteFile(path=str(i_new))

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        pattern = "*vcf*"
        dir_contents = list(dir.rglob(pattern))
        logger.info(
            f"Found following VCF files in {dir} matching {pattern}: {dir_contents}"
        )
        for fp in dir_contents:
            stem = str(fp.stem).split(".")[0]
            sample, caller = stem.split("_")[0:2]

            if sample not in samples:
                samples[sample] = VCF(sample=sample, caller=caller)

            if "tbi" in fp.name:
                setattr(samples[sample], "vcf_idx", FlyteFile(path=str(fp)))
            if "vcf" in fp.name and "tbi" not in fp.name:
                setattr(samples[sample], "vcf", FlyteFile(path=str(fp)))

        logger.info(f"Created following VCF objects from {dir}: {samples}")
        return list(samples.values())

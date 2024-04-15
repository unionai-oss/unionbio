from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from pathlib import Path


@dataclass
class VCF(DataClassJSONMixin):
    """
    Represents a VCF (Variant Call Format) file and its associated sample and index.

    This class defines the structure for representing a SAM file along with attributes
    that describe the associated sample and report.

    Attributes:
        sample (str): The name or identifier of the sample to which the SAM file belongs.
        aligner (str): The name of the aligner used to generate the SAM file.
        alignment (FlyteFile): A FlyteFile object representing the path to the alignment file.
        alignment_report (FlyteFile): A FlyteFile object representing an associated report
            for performance of the aligner.
        sorted (bool): A boolean value indicating whether the SAM file has been sorted.
        deduped (bool): A boolean value indicating whether the SAM file has been deduplicated.
        bqsr_report (FlyteFile): A FlyteFile object representing a report from the Base Quality
            Score Recalibration (BQSR) process.
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
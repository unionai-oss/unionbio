from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from pathlib import Path


@dataclass
class Alignment(DataClassJSONMixin):
    """
    Represents a SAM (Sequence Alignment/Map) file and its associated sample and report.

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
    aligner: str
    format: str
    alignment: FlyteFile | None = None
    alignment_idx: FlyteFile | None = None
    alignment_report: FlyteFile | None = None
    sorted: bool | None = None
    deduped: bool | None = None
    bqsr_report: FlyteFile | None = None

    def _get_state_str(self):
        state = f"{self.sample}_{self.aligner}"
        if self.sorted:
            state += "_sorted"
        if self.deduped:
            state += "_deduped"
        return state

    def get_alignment_fname(self):
        return f"{self._get_state_str()}_aligned.{self.format}"

    def get_alignment_idx_fname(self):
        return f"{self._get_state_str()}_aligned.bam.BAI"

    def get_report_fname(self):
        return f"{self._get_state_str()}_aligned_report.txt"

    def get_bqsr_fname(self):
        return f"{self._get_state_str()}_bqsr.table"

    def get_metrics_fname(self):
        return f"{self._get_state_str()}_metrics.txt"

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        for fp in list(dir.rglob("*aligned*")):
            sample, aligner = fp.stem.split("_")[0:2]

            if sample not in samples:
                samples[sample] = Alignment(sample=sample, aligner=aligner)

            if "sorted" in fp.name:
                setattr(samples[sample], "sorted", True)
            else:
                setattr(samples[sample], "sorted", False)

            if "deduped" in fp.name:
                setattr(samples[sample], "deduped", True)
            else:
                setattr(samples[sample], "deduped", False)

            if "bam" in fp.name or "sam" in fp.name:
                setattr(samples[sample], "alignment", FlyteFile(path=str(fp)))
            elif "report" in fp.name:
                setattr(samples[sample], "report", FlyteFile(path=str(fp)))

        return list(samples.values())

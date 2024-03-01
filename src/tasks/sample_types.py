from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from typing import Optional
from flytekit.types.file import FlyteFile
from config import logger
from pathlib import Path


@dataclass
class RawSample(DataClassJSONMixin):
    """
    Represents a raw sequencing sample via its associated files.

    This class defines the structure for representing a raw sequencing sample. It includes
    attributes for the sample name and paths to the raw read files (R1 and R2).

    Attributes:
        sample (str): The name or identifier of the raw sequencing sample.
        raw_r1 (FlyteFile): A FlyteFile object representing the path to the raw R1 read file.
        raw_r2 (FlyteFile): A FlyteFile object representing the path to the raw R2 read file.
    """

    sample: str
    raw_r1: FlyteFile | None = None
    raw_r2: FlyteFile | None = None

    def make_filenames(self):
        # Make filenames for filtered reads and report
        return (
            f"{self.sample}_1.fastq.gz",
            f"{self.sample}_2.fastq.gz",
        )

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        for fp in list(dir.rglob("*.fastq.gz")):
            logger.debug(f"Processing {fp}")
            sample, mate = fp.stem.strip(".fastq.gz").split("_")[0:2]
            logger.debug(f"Found sample {sample} and mate {mate}")
            if sample not in samples:
                samples[sample] = RawSample(sample=sample)
            if mate == "1":
                setattr(samples[sample], "raw_r1", FlyteFile(path=str(fp)))
            elif mate == "2":
                setattr(samples[sample], "raw_r2", FlyteFile(path=str(fp)))
        logger.info(f"Created {samples} from {dir}")
        return list(samples.values())


@dataclass
class FiltSample(DataClassJSONMixin):
    """
    Represents a filtered sequencing sample with its associated files and a quality report.

    This class defines the structure for representing a filtered sequencing sample. It includes
    attributes for the sample name, paths to the filtered read files (R1 and R2), and a quality
    report.

    Attributes:
        sample (str): The name or identifier of the filtered sequencing sample.
        filt_r1 (FlyteFile): A FlyteFile object representing the path to the filtered R1 read file.
        filt_r2 (FlyteFile): A FlyteFile object representing the path to the filtered R2 read file.
        report (FlyteFile): A FlyteFile object representing the quality report associated with
            the filtered sample.
    """

    sample: str
    filt_r1: FlyteFile | None = None
    filt_r2: FlyteFile | None = None
    report: FlyteFile | None = None

    def make_filenames(self):
        # Make filenames for filtered reads and report
        return (
            f"{self.sample}_1.filt.fastq.gz",
            f"{self.sample}_2.filt.fastq.gz",
            f"{self.sample}_filt-report.json",
        )

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        for fp in list(dir.rglob("*filt*")):
            if "fastq.gz" in fp.name:
                sample, mate = fp.stem.strip(".filt.fastq.gz").split("_")[0:2]
                logger.debug(f"Found sample {sample} and mate {mate} for {fp}")
            elif "report" in fp.name:
                sample = fp.stem.split("_")[0]
                mate = 0
                logger.debug(f"Found sample {sample} for report at {fp}")

            if sample not in samples:
                samples[sample] = FiltSample(sample=sample)

            if mate == "1":
                setattr(samples[sample], "filt_r1", FlyteFile(path=str(fp)))
            elif mate == "2":
                setattr(samples[sample], "filt_r2", FlyteFile(path=str(fp)))
            elif mate == 0:
                setattr(samples[sample], "report", FlyteFile(path=str(fp)))

        return list(samples.values())


@dataclass
class Alignment(DataClassJSONMixin):
    """
    Represents a SAM (Sequence Alignment/Map) file and its associated sample and report.

    This class defines the structure for representing a SAM file along with attributes
    that describe the associated sample and report.

    Attributes:
        sample (str): The name or identifier of the sample to which the SAM file belongs.
        aligner (str): The name of the aligner used to generate the SAM file.
        sam (FlyteFile): A FlyteFile object representing the path to the SAM file.
        alignment_report (FlyteFile): A FlyteFile object representing an associated report
            for performance of the aligner.
        sorted (bool): A boolean value indicating whether the SAM file has been sorted.
        deduped (bool): A boolean value indicating whether the SAM file has been deduplicated.
        bqsr_report (FlyteFile): A FlyteFile object representing a report from the Base Quality
            Score Recalibration (BQSR) process.
    """

    sample: str
    aligner: str
    sam: FlyteFile | None = None
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
        return f"{self._get_state_str()}_aligned.sam"

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

            if "sam" in fp.name:
                setattr(samples[sample], "sam", FlyteFile(path=str(fp)))
            elif "report" in fp.name:
                setattr(samples[sample], "report", FlyteFile(path=str(fp)))

        return list(samples.values())

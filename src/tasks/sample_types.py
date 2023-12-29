from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory
from flytekit import task
from config import base_image, logger
from pathlib import Path
from typing import List, Union, Any


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
    raw_r1: FlyteFile = FlyteFile(path="/dev/null")
    raw_r2: FlyteFile = FlyteFile(path="/dev/null")

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
    filt_r1: FlyteFile = FlyteFile(path="/dev/null")
    filt_r2: FlyteFile = FlyteFile(path="/dev/null")
    report: FlyteFile = FlyteFile(path="/dev/null")

    def make_filenames(self):
        # Make filenames for filtered reads and report
        return (
            f"{self.sample}_1_filt.fastq.gz",
            f"{self.sample}_2_filt.fastq.gz",
            f"{self.sample}_report.json",
        )

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        for fp in list(dir.rglob("*filt*")):
            if "fastq.gz" in fp.name:
                sample, mate = fp.stem.strip("fastq.gz").split("_")[0:2]
            elif "report" in fp.name:
                sample = fp.stem.split("_")[0]
                mate = 0

            if sample not in samples:
                samples[sample] = FiltSample(sample=sample)

            if mate == "1":
                setattr(samples[sample], "filt_r1", FlyteFile(path=str(fp)))
            elif mate == "2":
                setattr(samples[sample], "filt_r2", FlyteFile(path=str(fp)))

        return list(samples.values())


@dataclass
class SamFile(DataClassJSONMixin):
    """
    Represents a SAM (Sequence Alignment/Map) file and its associated sample and report.

    This class defines the structure for representing a SAM file along with attributes
    that describe the associated sample and report.

    Attributes:
        sample (str): The name or identifier of the sample to which the SAM file belongs.
        sam (FlyteFile): A FlyteFile object representing the path to the SAM file.
        report (FlyteFile): A FlyteFile object representing an associated report
            for performance of the aligner.
    """

    sample: str
    aligner: str
    sam: FlyteFile = FlyteFile(path="/dev/null")
    report: FlyteFile = FlyteFile(path="/dev/null")

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        for fp in list(dir.rglob("*aligned*")):
            sample, aligner = fp.stem.split("_")[0:2]

            if sample not in samples:
                samples[sample] = SamFile(sample=sample, aligner=aligner)

            if "sam" in fp.name:
                setattr(samples[sample], "sam", FlyteFile(path=str(fp)))
            elif "report" in fp.name:
                setattr(samples[sample], "report", FlyteFile(path=str(fp)))

        return list(samples.values())

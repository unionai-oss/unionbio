from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from pathlib import Path
from unionbio.config import logger


@dataclass
class Reads(DataClassJSONMixin):
    """
    Represents a sequencing reads sample via its associated fastq files.

    This class defines the structure for representing a sequencing sample. It includes
    attributes for the sample name and paths to the read files (R1 and R2).

    Attributes:
        sample (str): The name or identifier of the raw sequencing sample.
        filtered (bool): A boolean value indicating whether the reads have been filtered.
        filt_report (FlyteFile): A FlyteFile object representing the path to the filter report.
        read1 (FlyteFile): A FlyteFile object representing the path to the raw R1 read file.
        read2 (FlyteFile): A FlyteFile object representing the path to the raw R2 read file.
    """

    sample: str
    filtered: bool | None = None
    filt_report: FlyteFile | None = None
    read1: FlyteFile | None = None
    read2: FlyteFile | None = None

    def get_read_fnames(self):
        filt = "filt." if self.filtered else ""
        return (
            f"{self.sample}_1.{filt}fastq.gz",
            f"{self.sample}_2.{filt}fastq.gz",
        )

    def get_report_fname(self):
        return f"{self.sample}_fastq-filter-report.json"
    

    @classmethod
    def make_all(cls, dir: Path):
        samples = {}
        for fp in list(dir.rglob("*fastq*")):
            logger.debug(f"Processing {fp}")
            sample = fp.name.split("_")[0]
            logger.debug(f"Found sample {sample}")

            if sample not in samples:
                samples[sample] = Reads(sample=sample)

            if ".fastq.gz" in fp.name:
                mate = fp.name.strip(".fastq.gz").strip(".filt").split("_")[-1]
                logger.debug(f"Found mate {mate} for {sample}")
                if "1" in mate:
                    setattr(samples[sample], "read1", FlyteFile(path=str(fp)))
                elif "2" in mate:
                    setattr(samples[sample], "read2", FlyteFile(path=str(fp)))
            elif "filter-report" in fp.name:
                logger.debug(f"Found filter report for {sample}")
                setattr(samples[sample], "filtered", True)
                setattr(samples[sample], "filt_report", FlyteFile(path=str(fp)))

        logger.info(f"Created {samples} from {dir}")
        return list(samples.values())

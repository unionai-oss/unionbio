from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile


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
    raw_r1: FlyteFile
    raw_r2: FlyteFile


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
    filt_r1: FlyteFile
    filt_r2: FlyteFile
    report: FlyteFile


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
    sam: FlyteFile
    report: FlyteFile

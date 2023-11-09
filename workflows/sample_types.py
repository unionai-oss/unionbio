from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile


@dataclass
class RawSample(DataClassJSONMixin):
    sample: str
    raw_r1: FlyteFile
    raw_r2: FlyteFile


@dataclass
class FiltSample(DataClassJSONMixin):
    sample: str
    filt_r1: FlyteFile
    filt_r2: FlyteFile
    report: FlyteFile


@dataclass
class SamFile(DataClassJSONMixin):
    sample: str
    sam: FlyteFile
    report: FlyteFile

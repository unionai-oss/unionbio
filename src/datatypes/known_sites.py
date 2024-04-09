from pathlib import Path
from mashumaro.mixins.json import DataClassJSONMixin
from dataclasses import dataclass
from flytekit.types.file import FlyteFile
from config import logger


@dataclass
class Sites(DataClassJSONMixin):
    """
    Represents a known sites file and associated index file.

    Attributes:
        sites (FlyteFile): A FlyteFile object representing the path to the known sites file.
        idx (FlyteFile): A FlyteFile object representing the path to the known sites index file.
    """

    sites: FlyteFile
    idx: FlyteFile
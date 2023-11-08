import os
import subprocess
import logging
import gzip
import hashlib
import shutil
from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from typing import List, Tuple
from dataclasses import dataclass, asdict
from dataclasses_json import dataclass_json
from pathlib import Path
from mashumaro.mixins.json import DataClassJSONMixin

cache_repro = ShellTask(
    name="cache_repro",
    debug=True,
    cache=True,
    cache_version="1.0",
    script=
    """
    echo "Hello World, {inputs.msg}" > {outputs.o}
    """,
    inputs=kwtypes(msg=str),
    output_locs=[OutputLocation(var="o", var_type=FlyteFile, location='/root/repro.txt')],
)

@task(cache=True, cache_version="1.0")
def py_cache_repro(msg: str) -> FlyteFile:
    with open('dis.txt', 'w') as f:
        f.write(msg)
    return FlyteFile(path='dis.txt')
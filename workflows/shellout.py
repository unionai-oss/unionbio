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

@task#(container_image='localhost:30000/variant-discovery:latest')
def shelly():
# Create a subprocess to run the other Python script
    try:
        # Run the subprocess and capture the standard error
        result = subprocess.run(['hisat2'], check=True, text=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        # Format the traceback
        raise Exception(f"Command: {e.cmd}\nFailed with return code {e.returncode}:\n{e.stderr}")
    
    except FileNotFoundError as e:
        raise Exception(f"Process failed because the executable could not be found. Did you specify a container image in the task definition if using custom dependencies?\n{e}")
    
    else:
        # If there were no errors, print the output of the subprocess
        print("Subprocess finished successfully")
        print(result.stdout)
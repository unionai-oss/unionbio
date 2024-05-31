import os
import re
import subprocess as sp
from pathlib import Path
from flytekit import ImageSpec, task
from unionbio.config import folding_img
from flytekit.image_spec.image_spec import ImageBuildEngine

# Build package
result = sp.run(['poetry', 'build'], capture_output=True, text=True)

# Check if the command was successful
if result.returncode != 0:
    print(f"Poetry build failed with exit code {result.returncode}")
    print(result.stderr)
    exit(1)

# Extract the wheel name from the output
output = result.stdout
matches = re.findall(r'unionbio.*\.whl', output)
if matches and len(matches) == 1:
    wheel = matches[0]
else:
    print("Wheel file name not found in the output")
    exit(1)

folding_img.with_packages(f"dist/{wheel}")

ImageBuildEngine().build(folding_img)

# Parse output to grab image tag and write to config file somewhere
# Write task container_image to use the image tag
# Write poetry tasks to run tests using latest tag
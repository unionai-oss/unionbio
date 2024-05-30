import os
import subprocess as sp
from pathlib import Path
from flytekit import ImageSpec, task
from unionbio.config import test_spec, folding_img
from flytekit.image_spec.image_spec import ImageBuildEngine

# Build package
result = sp.run(["poetry", "build"])

folding_img.with_packages("dist/unionbio-0.1.0-py3-none-any.whl")
# ^^ Works

# ImageBuildEngine().build(folding_img)
# ^^ doesn't work, error: Distribution not found at: file:///dist/unionbio-0.1.0-py3-none-any.whl

@task(container_image=folding_img)
def foo():
    ...
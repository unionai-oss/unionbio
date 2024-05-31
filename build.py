import os
import re
import subprocess as sp
from pathlib import Path
from flytekit import ImageSpec, task
from unionbio.config import main_img, folding_img, parabricks_img
from flytekit.image_spec.image_spec import ImageBuildEngine

def built_wheel(output_dirname: str='dist', fmt: str='wheel') -> str:
    """Build the wheel package using Poetry and return the wheel file name."""

    # Build package
    result = sp.run(['poetry', 'build', '-o', output_dirname, '-f', fmt], capture_output=True, text=True)

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
    
    return f"{output_dirname}/{wheel}"

def update_img_config(cfg: Path, tags: dict[str, str]):
    
    with open(cfg, 'r') as f:
        cfg_content = f.read()
    
    for var, tag in tags.items():
        cfg_content = re.sub(rf"{var} = .+", f"{var} = \"{tag}\"", cfg_content)

    with open(cfg, 'w') as f:
        f.write(cfg_content)


# Run builds
whl = built_wheel()
config_path = Path('unionbio/config.py')
main_img.with_packages(whl)
folding_img.with_packages(whl)
parabricks_img.with_packages(whl)
tags = {
    'main_img_fqn': main_img.image_name(),
    'folding_img_fqn': folding_img.image_name(),
    'parabricks_img_fqn': parabricks_img.image_name()
}
update_img_config(config_path, tags)
ImageBuildEngine().build(main_img)
ImageBuildEngine().build(folding_img)
ImageBuildEngine().build(parabricks_img)

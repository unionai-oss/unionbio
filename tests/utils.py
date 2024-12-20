import os
import shutil
import filecmp
from pathlib import Path
from unionbio.config import logger

def comp_lines(file1, file2):
    with open(file1, "r") as f1, open(file2, "r") as f2:
        for line1, line2 in zip(f1, f2):
            if line1 != line2:
                logger.error("Lines do not match:")
                logger.error(line1.rstrip())
                logger.error(line2.rstrip())
                return False
        return True

def compare_dirs(expected: Path, actual: Path, mode: str) -> bool:
    
    assert mode in ["exists", "whole", "lines"], f"Please choose a valid comparison mode: exists, whole, lines"

    logger.debug(f"Comparing expected dir {expected} with test dir {actual} in mode {mode}")
    exp_files = [f for f in expected.iterdir() if f.is_file()]
    exp_dirs = [d for d in expected.iterdir() if d.is_dir()]

    for exp_file in exp_files:
        ac_file = actual.joinpath(exp_file.name)
        if mode == "exists":
            assert ac_file.exists(), f"Expected file {exp_file} does not exist in test directory"
        elif mode == "whole":
            assert filecmp.cmp(
                exp_file, ac_file
            ), f"{exp_file} and {ac_file} do not match"
        elif mode == "lines":
            assert comp_lines(exp_file, ac_file), f"Expected file {exp_file} and test file {ac_file} failed line by line comparison"

    for exp_dir in exp_dirs:
        ac_dir = actual.joinpath(exp_dir.name)
        assert ac_dir.exists(), f"Expected subdirectory {ac_dir} does not exist in test directory"
        compare_dirs(ac_dir, exp_dir, mode=mode)

    return True

def copy_dir_conts(src_dir, dest_dir):
    # Ensure the destination directory exists
    os.makedirs(dest_dir, exist_ok=True)

    # List all files and directories in the source directory
    items = os.listdir(src_dir)

    for item in items:
        src_path = os.path.join(src_dir, item)
        dest_path = os.path.join(dest_dir, item)

        # If the item is a directory, use copytree
        if os.path.isdir(src_path):
            shutil.copytree(src_path, dest_path, dirs_exist_ok=True)
        else:
            shutil.copy2(src_path, dest_path)



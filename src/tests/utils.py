import os
import filecmp


def dir_contents_match(dir1, dir2):
    dcmp = filecmp.dircmp(dir1, dir2)

    common_files = dcmp.common_files
    common_dirs = dcmp.common_dirs

    for common_file in common_files:
        file1_path = os.path.join(dir1, common_file)
        file2_path = os.path.join(dir2, common_file)
        assert filecmp.cmp(
            file1_path, file2_path
        ), f"{file1_path} and {file2_path} do not match"

    for common_dir in common_dirs:
        subdir1 = os.path.join(dir1, common_dir)
        subdir2 = os.path.join(dir2, common_dir)
        dir_contents_match(subdir1, subdir2)

    return True

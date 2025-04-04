import os
import subprocess
from pathlib import Path
from unionbio.types.protein import Protein
from unionbio.tasks.colabfold import cf_search, af_predict_local
from tests.utils import copy_dir_conts, compare_dirs
from tests.config import test_assets


# Check for GPU
def check_gpu():
    try:
        subprocess.check_output("nvidia-smi")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


assert check_gpu()


def test_cf_search(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    prot_in = Protein.make_all(tmp_path.joinpath("cf_search", "sequences"))[0]
    db_path = str(tmp_path.joinpath("cf_search", "databases"))
    prot_out = cf_search(
        prot=prot_in,
        db_path=db_path,
        search_args=[
            "--db1",
            "minisprot",
            "--use-env",
            "0",
            "--use-templates",
            "1",
            "--db2",
            "minisprot",
        ],
    )
    assert isinstance(prot_out, Protein)
    assert Path(prot_out.msa.path).name in os.listdir(
        tmp_path.joinpath("cf_search", "search_out")
    )


def test_af_predict_local(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    prot_in = Protein.make_all(tmp_path.joinpath("af_predict"))[0]
    predict_out = af_predict_local(
        prot=prot_in,
        mmcif_loc=str(tmp_path.joinpath("af_predict", "cifs")),
        batch_args=[
            "--num-models",
            "1",
            "--num-recycle",
            "1",
            "--stop-at-score",
            "1",
        ],
    )
    assert isinstance(predict_out, Protein)
    assert compare_dirs(
        tmp_path.joinpath("af_predict", "predict_out"),
        Path(predict_out.predict_out.path),
        mode="exists",
    )

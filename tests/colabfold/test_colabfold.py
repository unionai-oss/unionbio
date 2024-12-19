import os
from pathlib import Path
from unionbio.types.protein import Protein
from unionbio.tasks.colabfold import cf_search, af_predict
from tests.utils import copy_dir_conts
from tests.config import test_assets


def test_cf_search(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    prot_in = Protein.make_all(tmp_path.joinpath("cf_search", "sequences"))[0]
    db_path = str(tmp_path.joinpath("cf_search", "databases"))
    prot_out = cf_search(prot=prot_in, db_path=db_path, search_args=[
        "--db1", 
        "minisprot", 
        "--use-env",
        "0",
        "--use-templates",
        "1",
        "--db2",
        "minisprot",
        ]
    )
    assert isinstance(prot_out, Protein)
    assert Path(prot_out.msa.path).name in os.listdir(tmp_path.joinpath("cf_search", "search_out"))


def test_af_predict(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    prot_in = Protein.make_all(tmp_path.joinpath("af_predict"))[0]
    predict_out = af_predict(prot=prot_in, mmcif_loc=tmp_path.joinpath("af_predict", "cifs"), af_args=[
        "--templates",
        "--pdb-hit-file",
        "prot.hitfile.path",
        "--num-models",
        "1",
        "--num-recycle",
        "1",
        "--stop-at-score",
        "1",
    ])
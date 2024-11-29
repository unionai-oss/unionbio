from unionbio.types.protein import Protein
from unionbio.tasks.colabfold import cf_search, af_predict
from tests.utils import copy_dir_conts
from tests.config import test_assets


def test_cf_search(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    prot_in = Protein.make_all(tmp_path)
    db_path = str(tmp_path.joinpath("databases"))
    prot_out = cf_search(prot=prot, db_path=db_path, search_args=["--db1", "minisprot", "--use-env", "0"])
    # assert msa.path in os.listdir(tmp_path.joinpath("msas"))
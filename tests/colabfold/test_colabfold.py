# import sys
# sys.path.insert(0, '/home/flytekit/workspace/unionbio')

from flytekit.types.file import FlyteFile
from unionbio.tasks.colabfold import cf_search, af_predict
from tests.utils import copy_dir_conts
from tests.config import test_assets


def test_cf_search(tmp_path):
    copy_dir_conts(test_assets["protein_path"], tmp_path)
    seq = FlyteFile(path=str(tmp_path.joinpath("sequences", "frog.fasta")))
    db_path = str(tmp_path.joinpath("databases"))
    hf, msa = cf_search(seq=seq, db_path=db_path, search_args=["--db1", "minisprot", "--use-env", "0"])
    assert msa.path in os.listdir(tmp_path.joinpath("msas"))
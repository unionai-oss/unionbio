from unionbio.datatypes.reads import Reads
from unionbio.datatypes.protein import Protein
from pathlib import Path
from tasks.folding import prodigal_predict
from tests.config import test_assets

def test_prodigal_predict():
    reads = Reads.make_all(Path(test_assets["seq_dir"]).joinpath("folding"))[0]
    prot = prodigal_predict(reads)
    assert isinstance(prot, Protein)
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.protein import Protein
from pathlib import Path
from filecmp import cmp
from unionbio.tasks.folding import prodigal_predict
from tests.config import test_assets

def test_prodigal_predict():
    reads = Reads.make_all(Path(test_assets["folding_seq_dir"]))[0]
    prot = prodigal_predict(in_seq=reads)
    assert isinstance(prot, Protein)
    assert cmp(test_assets["prot_path"], prot.protein.path)
    assert cmp(test_assets["genes_path"], prot.genes.path)


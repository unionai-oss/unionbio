from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.utils import prepare_samples

def test_prep_samples_wf():
    samps = prepare_samples(seq_dir=FlyteDirectory(path="/Users/pryceturner/Desktop/genomic_data/sequence"))
    print(samps)
    assert samps == []
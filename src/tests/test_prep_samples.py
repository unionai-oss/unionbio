from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.sample_types import RawSample, prepare_raw_samples

def test_prep_samples_wf():
    samps = prepare_raw_samples(seq_dir=FlyteDirectory(path="/Users/pryceturner/Desktop/genomic_data/sequence"), sample_type='raw')
    print(samps)
    assert len(samps) == 2
    for i in samps:
        print(i.sample)
        assert isinstance(i, RawSample)
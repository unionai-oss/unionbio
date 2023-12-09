from pathlib import Path
from flytekit import workflow
from flytekit.types.directory import FlyteDirectory
from tasks.utils import prepare_raw_samples
from tasks.sample_types import RawSample, FiltSample, SamFile

def test_raw_sample_make_all():
    samps = RawSample.make_all(Path("/Users/pryceturner/Desktop/genomic_data/sequence"))
    assert len(samps) == 2
    for i in samps:
        assert isinstance(i, RawSample)

def test_filt_sample_fname():
    o1, o2, rep = FiltSample("test").make_filenames()
    assert o1 == "test_1_filt.fastq.gz"

def test_filt_sample_make_all():
    filt_samps = FiltSample.make_all(Path("/Users/pryceturner/Desktop/genomic_data/filtered_sequences"))
    assert len(filt_samps) == 2
    for i in filt_samps:
        assert isinstance(i, FiltSample)

def test_sam_file_make_all():
    sams = SamFile.make_all(Path("/Users/pryceturner/Desktop/genomic_data/alignments"))
    assert len(sams) == 2
    for i in sams:
        assert isinstance(i, SamFile)
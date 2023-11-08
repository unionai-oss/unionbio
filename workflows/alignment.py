import logging
from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task, dynamic
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.experimental import map_task
from typing import List, Tuple
from pathlib import Path

from .config import ref_loc
from .sample_types import FiltSample, SamFile
from .fastqc import fastqc
from .fastp import pyfastp
from .utils import prepare_samples
from .bowtie2 import bowtie2_align_paired_reads, bowtie2_index
from .hisat2 import hisat2_align_paired_reads, hisat2_index
from .multiqc import prep_multiqc_ins, multiqc, render_multiqc

@dynamic
def compare_aligners(bt2_idx: FlyteDirectory, hs2_idx: FlyteDirectory, samples: List[FiltSample]) -> List[List[SamFile]]:
    sams = []
    for sample in samples:
        bt2_sam = bowtie2_align_paired_reads(idx=bt2_idx, fs=sample)
        hs2_sam = hisat2_align_paired_reads(idx=hs2_idx, fs=sample)
        pair = [bt2_sam, hs2_sam]
        sams.append(pair)
    return sams

@workflow
def alignment_wf(seq_dir: FlyteDirectory='s3://my-s3-bucket/my-data/single'):
    qc = fastqc(seq_dir=seq_dir)
    samples = prepare_samples(seq_dir=seq_dir)
    filtered_samples = map_task(pyfastp)(rs=samples)
    bowtie2_idx = bowtie2_index(ref=ref_loc)
    hisat2_idx = hisat2_index(ref=ref_loc)
    sams = compare_aligners(bt2_idx=bowtie2_idx, hs2_idx=hisat2_idx, samples=filtered_samples)
    mqc_prep = prep_multiqc_ins(fqc=qc, filt_reps=filtered_samples, sams=sams)
    mqc = multiqc(report_dir=mqc_prep)
    render_multiqc(report=mqc)
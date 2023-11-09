import logging
from datetime import timedelta
from flytekit import kwtypes, workflow, ImageSpec, Resources, current_context, task, dynamic, approve
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.experimental import map_task
from typing import List, Tuple
from pathlib import Path

from .config import ref_loc, seq_dir
from .sample_types import FiltSample, SamFile
from .fastqc import fastqc
from .fastp import pyfastp
from .utils import prepare_samples, pass_qc
from .bowtie2 import bowtie2_align_paired_reads, bowtie2_index
from .hisat2 import hisat2_align_paired_reads, hisat2_index
from .multiqc import render_multiqc

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
def alignment_wf(seq_dir: FlyteDirectory=seq_dir) -> FlyteFile:
    fqc_dir = fastqc(seq_dir=seq_dir)
    samples = prepare_samples(seq_dir=seq_dir)
    filtered_samples = map_task(pyfastp)(rs=samples)
    
    qc_check = render_multiqc(fqc=fqc_dir, filt_reps=filtered_samples, sams=[])
    report = approve(qc_check, "qc-report-approval", timeout=timedelta(hours=2))

    bowtie2_idx = bowtie2_index(ref=ref_loc)
    hisat2_idx = hisat2_index(ref=ref_loc)
    sams = compare_aligners(bt2_idx=bowtie2_idx, hs2_idx=hisat2_idx, samples=filtered_samples)

    report >> bowtie2_idx
    report >> hisat2_idx
    report >> sams

    return render_multiqc(fqc=fqc_dir, filt_reps=filtered_samples, sams=sams)
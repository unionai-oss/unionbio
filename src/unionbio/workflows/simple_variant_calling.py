from flytekit import workflow, LaunchPlan
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit import map_task, task, dynamic
from typing import List
from unionbio.config import main_img_fqn, remote_reads, remote_ref, remote_sites_vcf, remote_sites_idx
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.datatypes.alignment import Alignment
from unionbio.tasks.fastqc import fastqc
from unionbio.tasks.fastp import pyfastp
from unionbio.tasks.utils import prepare_raw_samples, reformat_alignments, fetch_remote_reads, fetch_remote_reference, fetch_remote_sites
from unionbio.tasks.bowtie2 import bowtie2_align_samples, bowtie2_index
from unionbio.tasks.multiqc import render_multiqc
from unionbio.tasks.base_recal import recalibrate_samples
from unionbio.tasks.mark_dups import mark_dups_samples
from unionbio.tasks.sort_sam import sort_samples
from unionbio.tasks.haplotype_caller import hc_call_samples

@workflow
def calling_wf(
    ref_url: str = remote_ref,
    reads_urls: List[str] = remote_reads,
    remote_sites_vcf: str = remote_sites_vcf,
    remote_sites_idx: str = remote_sites_idx,
):# -> FlyteFile:
    """
    Run an alignment workflow on FastQ files contained in the configured seq_dir.

    This function performs QC, preprocessing, index generation, and finally alignments
    on FastQ input files present in a pre-configured S3 prefix. It returns a FlyteFile
    of a MultiQC report containing all relevant statistics of the different steps.

    Args:
        seq_dir (FlyteDirectory, optional): The input directory containing sequencing data.
            Defaults to the value of `seq_dir` (if provided).

    Returns:
        FlyteFile: A FlyteFile object representing the output of the alignment workflow.
    """
    # Fetch remote inputs
    ref = fetch_remote_reference(url=ref_url)
    reads = fetch_remote_reads(urls=reads_urls)
    sites = fetch_remote_sites(sites=remote_sites_vcf, idx=remote_sites_idx)

    # Generate FastQC reports and check for failures
    fqc_out = fastqc(reads=reads)

    # # Map out filtering across all samples and generate indices
    filtered_samples = map_task(pyfastp)(rs=reads)

    # # Explicitly define task dependencies
    fqc_out >> filtered_samples

    # # Generate a bowtie2 index or load it from cache
    bowtie2_idx = bowtie2_index(ref=ref)

    # # Generate alignments using bowtie2
    sams = bowtie2_align_samples(idx=bowtie2_idx, samples=filtered_samples)

    # # Recalibrate & Reformat
    # deduped = mark_dups_samples(als=sams)
    # sorted = sort_samples(als=deduped)
    # recal_sams = recalibrate_samples(als=sams, sites=sites, ref=ref)
    # bams = reformat_alignments(als=recal_sams, to_format='bam')

    # # Call Variants
    # vcfs = hc_call_samples(ref=ref_path, als=sams)

    # # Generate final multiqc report with stats from all steps
    # rep = render_multiqc(fqc=fqc_out, filt_reps=filtered_samples, sams=sams, vcfs=vcfs)

    # return rep

# call_lp = LaunchPlan.get_or_create(calling_wf, name="calling_wf_exp", default_inputs={"seq_dir": "s3://my-s3-bucket/my-data/sequences", "ref_path": "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"})
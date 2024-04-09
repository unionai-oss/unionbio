from typing import List
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile

from tasks.utils import fetch_remote_reads, fetch_remote_reference, fetch_remote_sites
from tasks.bwa import bwa_index
from tasks.parabricks import fq2bam


@task
def index_reference(ref: str) -> FlyteDirectory:
    return FlyteDirectory(path="/tmp")


@task
def get_known_sites(sites: str, idx: str) -> FlyteDirectory:
    return FlyteDirectory(path="/tmp")


@task
def pb_deepvar(bam_dir: FlyteDirectory, ref_dir: FlyteDirectory) -> FlyteDirectory:
    return FlyteDirectory(path="/tmp")


@task
def pb_haplocall(
    bam_dir: FlyteDirectory, recal: FlyteFile, ref_dir: FlyteDirectory
) -> FlyteDirectory:
    return FlyteDirectory(path="/tmp")


@task
def intersect_vars(vcf1: FlyteDirectory, vcf2: FlyteDirectory) -> FlyteFile:
    return FlyteFile(path="/tmp")


@workflow
def call_vars(
    reads: List[str] = [
        'wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseqx/wgs_pcr_free/30x/HG002.hiseqx.pcr-free.30x.R1.fastq.gz',
        'wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseqx/wgs_pcr_free/30x/HG002.hiseqx.pcr-free.30x.R2.fastq.gz'
        ],
    ref: str = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    sites: List[str] = [
        'wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
        'wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
        ],
) -> FlyteFile:
    read_obj = fetch_remote_reads(urls=reads)
    ref_obj = fetch_remote_reference(url=ref)
    sites_obj = fetch_remote_sites(sites=sites[0], idx=sites[1])
    ref_idx = bwa_index(ref=ref_obj)
    alignment = fq2bam(reads=read_obj, sites=sites_obj, ref=ref_obj)
    deepvar_dir = pb_deepvar(bam_dir=bam_dir, ref_dir=ref_dir)
    haplocall_dir = pb_haplocall(bam_dir=bam_dir, recal=recal, ref_dir=ref_dir)
    return intersect_vars(vcf1=deepvar_dir, vcf2=haplocall_dir)


# @workflow
# def comparison_wf() -> typing.Tuple[bool, bool, str, str, str]:
#     data = get_data(
#         url="https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"
#     )
    # ff1, s1 = dgx_pb_align(indir=data)
    # ff2, s2 = dgx_basic_align(indir=data)
    # ff3, s3 = demo_basic_align(indir=data)
    # c1 = compare_bams(in1=ff1, in2=ff2)
    # c2 = compare_bams(in1=ff1, in2=ff3)
    # return c1, c2, s1, s2, s3

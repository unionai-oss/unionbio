from typing import List
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile

from tasks.utils import fetch_remote_reads, fetch_remote_reference, fetch_remote_sites, intersect_vcfs
from tasks.bwa import bwa_index
from tasks.parabricks import fq2bam, pb_deepvar, pb_haplocall


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
    alignment = fq2bam(reads=read_obj, sites=sites_obj, ref=ref_idx)
    deepvar_vcf = pb_deepvar(bam=alignment, ref=ref)
    haplocall_vcf = pb_haplocall(bam=alignment, ref=ref)
    return intersect_vcfs(vcf1=deepvar_vcf, vcf2=haplocall_vcf)


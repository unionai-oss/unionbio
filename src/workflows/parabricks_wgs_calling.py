import typing
from flytekit import task, workflow
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile

from tasks.utils import get_data


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
    r1: str = "read1.fq.gz",
    r2: str = "read2.fq.gz",
    ref: str = "hg38.fasta",
    sites: str = "known_sites.vcf.gz",
    sites_idx: str = "known_sites.vcf.gz.tbi",
) -> FlyteFile:
    read_dir = get_data(r1=r1, r2=r2)
    ref_dir = index_reference(ref=ref)
    sites_dir = get_known_sites(sites=sites, idx=sites_idx)
    bam_dir, recal = pb_fq2bam(reads=read_dir, ref_dir=ref_dir, sites=sites_dir)
    deepvar_dir = pb_deepvar(bam_dir=bam_dir, ref_dir=ref_dir)
    haplocall_dir = pb_haplocall(bam_dir=bam_dir, recal=recal, ref_dir=ref_dir)
    return intersect_vars(vcf1=deepvar_dir, vcf2=haplocall_dir)


@workflow
def comparison_wf() -> typing.Tuple[bool, bool, str, str, str]:
    data = get_data(
        url="https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"
    )
    # ff1, s1 = dgx_pb_align(indir=data)
    # ff2, s2 = dgx_basic_align(indir=data)
    # ff3, s3 = demo_basic_align(indir=data)
    # c1 = compare_bams(in1=ff1, in2=ff2)
    # c2 = compare_bams(in1=ff1, in2=ff3)
    # return c1, c2, s1, s2, s3

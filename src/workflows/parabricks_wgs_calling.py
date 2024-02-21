import os
import time
import typing
import requests
from pathlib import Path
from flytekit import task, workflow, ImageSpec, Resources
from flytekit.extras.tasks.shell import subproc_execute
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
# from flytekitplugins.dgx import DGXConfig

@task
def get_data(r1: str, r2: str) -> FlyteDirectory:
    return FlyteDirectory(path='/tmp')
    # r1p = Path(r1.split('/')[-1])
    # r2p = Path(r2.split('/')[-1])

@task
def index_reference(ref: str) -> FlyteDirectory:
    return FlyteDirectory(path='/tmp')

@task
def get_known_sites(sites: str, idx: str) -> FlyteDirectory:
    return FlyteDirectory(path='/tmp')

@task#(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=imageSpecNvcr)
def pb_fq2bam(reads: FlyteDirectory, ref_dir: FlyteDirectory, sites: FlyteDirectory) -> typing.Tuple[FlyteDirectory, FlyteFile]:
    
    outpath = 'out.bam'
    # indir.download()
    # loc_dir = Path(indir.path)
    # r1 = loc_dir.joinpath('Data/sample_1.fq.gz')
    # r2 = loc_dir.joinpath('Data/sample_2.fq.gz')
    # ref = loc_dir.joinpath('Ref/Homo_sapiens_assembly38.fasta')
    # sites = loc_dir.joinpath('Ref/Homo_sapiens_assembly38.known_indels.vcf.gz')
    recal_out = 'recal_data.table'

    # t1 = time.time()
    # out, err = subproc_execute([
    #     'pbrun',
    #     'fq2bam',
    #     '--ref',
    #     str(ref),
    #     '--in-fq',
    #     str(r1),
    #     str(r2),
    #     '--knownSites',
    #     str(sites),
    #     '--out-bam',
    #     outpath,
    #     '--out-recal-file',
    #     recal_out
    #     ])
    # elapsed = f'pbrun fq2bam in DGX took {time.time() - t1} seconds'

    return FlyteDirectory(path='/tmp'), FlyteFile(path='/dev/null')

@task
def pb_deepvar(bam_dir: FlyteDirectory, ref_dir: FlyteDirectory) -> FlyteDirectory:
    return FlyteDirectory(path='/tmp')

@task
def pb_haplocall(bam_dir: FlyteDirectory, recal: FlyteFile, ref_dir:FlyteDirectory) -> FlyteDirectory:
    return FlyteDirectory(path='/tmp')

@task
def intersect_vars(vcf1: FlyteDirectory, vcf2: FlyteDirectory) -> FlyteFile:
    return FlyteFile(path='/tmp')

@workflow
def call_vars(r1: str='read1.fq.gz', r2: str='read2.fq.gz', ref: str='hg38.fasta', sites: str='known_sites.vcf.gz', sites_idx: str='known_sites.vcf.gz.tbi') -> FlyteFile:
    read_dir = get_data(r1=r1, r2=r2)
    ref_dir = index_reference(ref=ref)
    sites_dir = get_known_sites(sites=sites, idx=sites_idx)
    bam_dir, recal = pb_fq2bam(reads=read_dir, ref_dir=ref_dir, sites=sites_dir)
    deepvar_dir = pb_deepvar(bam_dir=bam_dir, ref_dir=ref_dir)
    haplocall_dir = pb_haplocall(bam_dir=bam_dir, recal=recal, ref_dir=ref_dir)
    return intersect_vars(vcf1=deepvar_dir, vcf2=haplocall_dir)
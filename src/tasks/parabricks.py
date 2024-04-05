from typing import List, Tuple
from pathlib import Path
from flytekit import task, Resources, current_context
from flytekit.extras.tasks.shell import subproc_execute
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekitplugins.dgx import DGXConfig

from config import pb_image
from datatypes.alignment import Alignment


# Supported DGX instances:
# dgxa100.80g.1.norm, dgxa100.80g.2.norm, dgxa100.80g.4.norm, dgxa100.80g.8.norm
@task(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=pb_image)
def fq2bam(reads: List[FlyteFile], sites: List[FlyteFile], ref_name: str, ref_dir: FlyteDirectory) -> Tuple[FlyteFile, FlyteFile]:
    """
    Takes an input directory containing sequence data and an indexed reference genome and
    performs alignment using Parabricks' fq2bam tool.

    Args:
        

    Returns:
        
    """

    RGTAG = "@RG\tID:HG002\tLB:lib\tPL:Illumina\tSM:HG002\tPU:HG002"
    reads[0].download()
    reads[1].download()
    sites[0].download()
    sites[1].download()
    ref_dir.download()

    bam_out = "out.bam"
    recal_out = "recal_data.table"

    out, err = subproc_execute(
        [
            "pbrun",
            "fq2bam",
            "--ref",
            str(ref_name),
            "--in-fq",
            str(reads[0].path),
            str(reads[1].path),
            RGTAG,
            "--knownSites",
            str(sites),
            "--out-bam",
            bam_out,
            "--out-recal-file",
            recal_out,
        ]
    )

    

    return FlyteFile(path=bam_out), FlyteFile(path=recal_out)


@task(requests=Resources(gpu="1", mem="32Gi", cpu="32"), container_image=pb_image)
# @task(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=pb_image)
def basic_align(indir: FlyteDirectory) -> Tuple[FlyteFile, str]:
    """
    Aligns paired-end sequencing reads using BWA-MEM and GATK tools, and returns the path to the processed BAM file
    and the elapsed time for the alignment process.

    Args:
        indir (FlyteDirectory): The input directory containing sequencing read files and reference data.

    Returns:
        FlyteFile: The path to the processed BAM file in FlyteFile format.
    """
    indir.download()
    loc_dir = Path(indir.path)
    r1 = loc_dir.joinpath("Data/sample_1.fq.gz")
    r2 = loc_dir.joinpath("Data/sample_2.fq.gz")
    ref = loc_dir.joinpath("Ref/Homo_sapiens_assembly38.fasta")
    sites = loc_dir.joinpath("Ref/Homo_sapiens_assembly38.known_indels.vcf.gz")
    bampath = "bwa_mem_out.bam"

    cmd = " ".join(
        [
            "bwa",
            "mem",
            "-t",
            "32",
            "-K",
            "10000000",
            "-R",
            r'"@RG\tID:sample_rg1\tLB:lib1\tPL:bar\tSM:sample\tPU:sample_rg1"',
            str(ref),
            str(r1),
            str(r2),
            "|",
            "java",
            "-jar",
            "/usr/local/bin/gatk",
            "SortSam",
            "--MAX_RECORDS_IN_RAM",
            "5000000",
            "-I",
            "/dev/stdin",
            "-O",
            bampath,
            "--SORT_ORDER",
            "coordinate",
        ]
    )

    out1, err1 = subproc_execute(cmd, shell=True)

    dup_bam = "mark_dups.bam"
    out2, err2 = subproc_execute(
        [
            "java",
            "-jar",
            "/usr/local/bin/gatk",
            "MarkDuplicates",
            "-I",
            bampath,
            "-O",
            dup_bam,
            "-M",
            "metrics.txt",
        ]
    )

    recal_out = "recal_data.table"
    out3, err3 = subproc_execute(
        [
            "java",
            "-jar",
            "/usr/local/bin/gatk",
            "BaseRecalibrator",
            "--input",
            dup_bam,
            "--output",
            recal_out,
            "--known-sites",
            str(sites),
            "--reference",
            str(ref),
        ]
    )

    return FlyteFile(path=dup_bam)

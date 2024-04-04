import time
import typing
from pathlib import Path
from flytekit import task, Resources
from flytekit.extras.tasks.shell import subproc_execute
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekitplugins.dgx import DGXConfig

from config import pb_image


# Supported DGX instances:
# dgxa100.80g.1.norm, dgxa100.80g.2.norm, dgxa100.80g.4.norm, dgxa100.80g.8.norm
@task(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=pb_image)
def dgx_pb_align(indir: FlyteDirectory) -> typing.Tuple[FlyteFile, str]:
    """
    Takes an input directory containing sequence data and an indexed reference genome and
    performs alignment using Parabricks' fq2bam tool.

    Args:
        indir (FlyteDirectory): The input directory containing sequences in the Data sub-directory
        and an indexed reference genome in the Ref sub-directory.

    Returns:
        typing.Tuple[FlyteFile, str]: A tuple containing:
            - The aligned BAM file in FlyteFile format.
            - A string representing run statistics about the alignment process.
    """

    outpath = "out.bam"
    indir.download()
    loc_dir = Path(indir.path)
    r1 = loc_dir.joinpath("Data/sample_1.fq.gz")
    r2 = loc_dir.joinpath("Data/sample_2.fq.gz")
    ref = loc_dir.joinpath("Ref/Homo_sapiens_assembly38.fasta")
    sites = loc_dir.joinpath("Ref/Homo_sapiens_assembly38.known_indels.vcf.gz")
    recal_out = "recal_data.table"

    t1 = time.time()
    out, err = subproc_execute(
        [
            "pbrun",
            "fq2bam",
            "--ref",
            str(ref),
            "--in-fq",
            str(r1),
            str(r2),
            "--knownSites",
            str(sites),
            "--out-bam",
            outpath,
            "--out-recal-file",
            recal_out,
        ]
    )
    elapsed = f"pbrun fq2bam in DGX took {time.time() - t1} seconds"

    return FlyteFile(path=outpath), elapsed


@task(requests=Resources(gpu="1", mem="32Gi", cpu="32"), container_image=pb_image)
# @task(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=pb_image)
def basic_align(indir: FlyteDirectory) -> typing.Tuple[FlyteFile, str]:
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

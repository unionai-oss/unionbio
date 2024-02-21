import os
import time
import typing
import requests
import tarfile
from pathlib import Path
from flytekit import task, workflow, Resources, current_context
from flytekit.extras.tasks.shell import subproc_execute
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekitplugins.dgx import DGXConfig

pb_img = "ghcr.io/unionai/dgx-parabricks:20240202"


@task(container_image=pb_img)
def get_data(url: str) -> FlyteDirectory:
    """
    Downloads a tar.gz file from the specified URL, extracts its contents, and returns a FlyteDirectory object.

    Args:
        url (str): The URL of the tar.gz file to download.

    Returns:
        FlyteDirectory: A FlyteDirectory object representing the directory where the contents of the tar.gz file were extracted.

    Raises:
        requests.HTTPError: If an HTTP error occurs while downloading the file.
    """
    try:
        response = requests.get(url)
        tar_name = url.split("/")[-1]
        with open(tar_name, "wb") as file:
            file.write(response.content)
        working_dir = current_context().working_directory
    except requests.HTTPError as e:
        print(f"HTTP error: {e}")
        raise e

    out_dir = Path(os.path.join(working_dir, tar_name))
    with tarfile.open(tar_name, "r:gz") as tarf:
        tarf.extractall(out_dir)

    return FlyteDirectory(path=out_dir)


# Supported DGX instances:
# dgxa100.80g.1.norm, dgxa100.80g.2.norm, dgxa100.80g.4.norm, dgxa100.80g.8.norm
@task(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=pb_img)
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


@task(requests=Resources(gpu="1", mem="32Gi", cpu="32"), container_image=pb_img)
def demo_basic_align(indir: FlyteDirectory) -> typing.Tuple[FlyteFile, str]:
    return basic_align(indir=indir, env="Demo with 32Gi mem and 32 CPU")


@task(task_config=DGXConfig(instance="dgxa100.80g.1.norm"), container_image=pb_img)
def dgx_basic_align(indir: FlyteDirectory) -> typing.Tuple[FlyteFile, str]:
    return basic_align(indir=indir, env="DGX A100 80G")


def basic_align(indir: FlyteDirectory, env: str) -> typing.Tuple[FlyteFile, str]:
    """
    Aligns paired-end sequencing reads using BWA-MEM and GATK tools, and returns the path to the processed BAM file
    and the elapsed time for the alignment process.

    Args:
        indir (FlyteDirectory): The input directory containing sequencing read files and reference data.
        env (str): A description of the environment in which the alignment is performed.

    Returns:
        typing.Tuple[FlyteFile, str]: A tuple containing:
            - The path to the processed BAM file in FlyteFile format.
            - A string indicating the elapsed time for the alignment process and the environment description.
    """
    indir.download()
    loc_dir = Path(indir.path)
    r1 = loc_dir.joinpath("Data/sample_1.fq.gz")
    r2 = loc_dir.joinpath("Data/sample_2.fq.gz")
    ref = loc_dir.joinpath("Ref/Homo_sapiens_assembly38.fasta")
    sites = loc_dir.joinpath("Ref/Homo_sapiens_assembly38.known_indels.vcf.gz")
    bampath = "bwa_mem_out.bam"

    t1 = time.time()
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
            "pb_imgRecalibrator",
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

    elapsed = f"bwa + gatk in {env} took {time.time() - t1} seconds"

    return FlyteFile(path=dup_bam), elapsed


@task(container_image=pb_img)
def compare_bams(in1: FlyteFile, in2: FlyteFile) -> bool:
    """
    Compares two BAM files and returns True if they are identical, False otherwise.

    Args:
        in1 (FlyteFile): The first input BAM file.
        in2 (FlyteFile): The second input BAM file.

    Returns:
        bool: True if the BAM files are identical, False otherwise.
    """
    in1.download()
    in2.download()

    cmp1 = [
        "bam",
        "diff",
        "--in1",
        in1.path,
        "--in2",
        in2.path,
        "--noCigar",
        "--isize",
        "--flag",
        "--mate",
        "--mapQual",
    ]

    out, err = subproc_execute(cmp1)

    no_out = out == "" and err == ""

    return no_out


@workflow
def comparison_wf() -> typing.Tuple[bool, bool, str, str, str]:
    data = get_data(
        url="https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz"
    )
    ff1, s1 = dgx_pb_align(indir=data)
    ff2, s2 = dgx_basic_align(indir=data)
    ff3, s3 = demo_basic_align(indir=data)
    c1 = compare_bams(in1=ff1, in2=ff2)
    c2 = compare_bams(in1=ff1, in2=ff3)
    return c1, c2, s1, s2, s3

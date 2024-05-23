import os
import zipfile
import requests
import tarfile
import gzip
import ftplib
from pathlib import Path
from typing import List
from flytekit import task, current_context
from flytekit.types.directory import FlyteDirectory
from flytekit.types.file import FlyteFile
from flytekit.configuration import Config
from flytekit.remote import FlyteRemote
from flytekit.extras.tasks.shell import subproc_execute

from unionbio.config import main_img, logger, parabricks_img
from unionbio.datatypes.reads import Reads
from unionbio.datatypes.reference import Reference
from unionbio.datatypes.variants import VCF
from unionbio.tasks.helpers import fetch_file


@task(container_image=main_img)
def prepare_raw_samples(seq_dir: FlyteDirectory) -> List[Reads]:
    """
    Prepare and process raw sequencing data to create a list of RawSample objects.

    This function processes raw sequencing data located in the specified input directory
    and prepares it to create a list of RawSample objects.

    Args:
        seq_dir (FlyteDirectory): The input directory containing raw sequencing data.

    Returns:
        List[Reads]: A list of Reads objects representing the processed sequencing data.
    """
    seq_dir.download()
    return Reads.make_all(Path(seq_dir))


@task(cache=True, cache_version="1.0")
def fetch_remote_reference(url: str) -> Reference:
    workdir = current_context().working_directory
    ref_path = fetch_file(url, workdir)
    return Reference(ref_name=str(ref_path.name), ref_dir=FlyteDirectory(path=workdir))


@task(cache=True, cache_version="1.0")
def fetch_remote_reads(urls: List[str]) -> Reads:
    """
    Fetches remote reads from a list of URLs and returns a list of Reads objects.

    Args:

    Returns:

    """
    workdir = current_context().working_directory
    for url in urls:
        fetch_file(url, workdir)
    return Reads.make_all(Path(workdir))[0]


@task(cache=True, cache_version="1.0")
def fetch_remote_sites(sites: str, idx: str) -> VCF:
    """
    Fetches remote known sites from a URL and returns a Sites object.

    Args:

    Returns:

    """
    workdir = current_context().working_directory
    sites_path = fetch_file(sites, workdir)
    idx_path = fetch_file(idx, workdir)
    return VCF(sites=FlyteFile(path=sites_path), idx=FlyteFile(path=idx_path))


@task
def fetch_files(urls: List[str]) -> List[FlyteFile]:
    outfiles = []
    workdir = current_context().working_directory
    for url in urls:
        lpath = fetch_file(url, workdir)
        outfiles.append(FlyteFile(path=lpath))
    return outfiles


@task
def fetch_tarfile(url: str) -> FlyteDirectory:
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
    except requests.HTTPError as e:
        print(f"HTTP error: {e}")
        raise e

    working_dir = current_context().working_directory
    out_dir = Path(os.path.join(working_dir, tar_name))
    with tarfile.open(tar_name, "r:gz") as tarf:
        tarf.extractall(out_dir)

    return FlyteDirectory(path=out_dir)


@task
def check_fastqc_reports(rep_dir: FlyteDirectory) -> str:
    """
    Check FastQC reports for errors.

    This function checks FastQC reports for errors and raises an exception if any are found.

    Args:
        rep_dir (FlyteDirectory): The input directory containing FastQC reports.
    """
    rep_dir.download()
    all_zips = list(Path(rep_dir.path).rglob("*fastqc.zip*"))

    for p in all_zips:
        logger.debug(f"Checking {p}")
        with zipfile.ZipFile(p, "r") as zip_file:
            logger.debug(f"{zip_file.filename}")
            logger.debug(f"Archive contains {zip_file.namelist()}")
            with zip_file.open(
                f"{Path(zip_file.filename).stem}/summary.txt"
            ) as summary:
                contents = summary.read().decode("utf-8")
                logger.debug(f"Contents of summary.txt: {contents}")
                if "FAIL" in contents:
                    return "FAIL"
                elif "WARN" in contents:
                    return "WARN"

    return "PASS"


@task(container_image=parabricks_img)
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

    result = subproc_execute(cmp1)

    no_out = result.out == "" and result.err == ""

    return no_out


@task(container_image=main_img)
def intersect_vcfs(vcf1: VCF, vcf2: VCF) -> VCF:
    """
    Takes the intersection of 2 VCF files and returns a new VCF file to increase
    calling sensitivity.

    Args:
        vcf1 (VCF): The first input VCF object.
        vcf2 (VCF): The second input VCF object.

    Returns:
        VCF: Intersected and zipped VCF object.
    """
    vcf1.dl_all()
    vcf2.dl_all()
    isec_out = VCF(sample=vcf1.sample, caller=f"{vcf1.caller}_{vcf2.caller}_isec")
    fname_out = isec_out.get_vcf_fname()

    cmd = " ".join([
        "bcftools",
        "isec",
        "-n=2",
        "-O",
        "z",
        "-w",
        "1",
        vcf1.vcf.path,
        vcf2.vcf.path,
        "-o",
        fname_out,
    ])

    subproc_execute(cmd, shell=True)
    setattr(isec_out, "vcf", FlyteFile(path=fname_out))

    return isec_out

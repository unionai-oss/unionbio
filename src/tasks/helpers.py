import gzip
import ftplib
import requests
from pathlib import Path
from flytekit.remote import FlyteRemote
from flytekit.configuration import Config

def gunzip_file(gzip_file: Path) -> Path:
    # Ensure the input file exists
    if not gzip_file.exists():
        raise FileNotFoundError(f"{gzip_file} not found")

    # Ensure the file has a .gz extension
    if not gzip_file.suffix == '.gz':
        raise ValueError("Input file is not a gzip file")

    # Define the output file path
    output_file = gzip_file.with_suffix("")

    # Open input and output files
    with gzip.open(gzip_file, 'rb') as f_in, open(output_file, 'wb') as f_out:
        # Decompress and write to output file
        f_out.write(f_in.read())

    return output_file

def get_remote(local=None, config_file=None):
    """
    Get remote configuration settings and return a remote object.

    This function retrieves remote configuration settings, including the local flag and
    a configuration file, and uses them to create and return a remote object.

    Args:
        local (bool, optional): A flag indicating whether to use local settings. If True,
            the function will use local settings; if False, it will use remote settings.
            Defaults to None, which implies the use of default settings.
        config_file (str, optional): The path to a custom configuration file. If provided,
            this file will be used for configuration settings. Defaults to None.

    Returns:
        Remote: A remote object configured with the specified settings.
    """
    return FlyteRemote(
        config=Config.auto(
            config_file=(
                None
                if local
                else config_file
                if config_file is not None
                else str(Path.home() / ".flyte" / "config-sandbox.yaml")
            )
        ),
        default_project="flytesnacks",
        default_domain="development",
    )

def fetch_file(url: str, local_dir: str) -> Path:
    """
    Downloads a file from the specified URL.

    Args:
        url (str): The URL of the tar.gz file to download.
        local_dir (Path): The directory where you would like this file saved.

    Returns:
        Path: The local path to the file.

    Raises:
        requests.HTTPError: If an HTTP error occurs while downloading the file.
    """
    url_parts = url.split("/")
    host = url_parts[2]
    prot = url_parts[0]
    fname = url_parts[-1]
    remote_dir = "/".join(url_parts[3:-1])
    local_path = Path(local_dir).joinpath(fname)

    if prot == "ftp:":  # FTP
        ftp = ftplib.FTP(host)
        ftp.login()
        ftp.cwd(remote_dir)
        with open(local_path, "wb") as file:
            ftp.retrbinary(f"RETR {fname}", file.write)
        ftp.quit()
    elif prot == "http:" or prot == "https:":  # HTTP
        try:
            response = requests.get(url)
            with open(local_path, "wb") as file:
                file.write(response.content)
        except requests.HTTPError as e:
            print(f"HTTP error: {e}")
            raise e
    return local_path
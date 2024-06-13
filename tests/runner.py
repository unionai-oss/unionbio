import docker
import os
from config import proj_rt, main_img_test_fqn

def run_pytest_in_docker():
    client = docker.from_env()
    
    # try:
    # Run the Docker container
    print("\nRunning pytest in Docker container...")
    container = client.containers.run(
        image=main_img_test_fqn,
        volumes={
            proj_rt.joinpath('tests'): {
                'bind': '/root/tests',
                'mode': 'rw',
            },
            # Mount unionbio so it's discoverable in PYTHONPATH
            proj_rt.joinpath("src/unionbio"): {
                'bind': '/root/unionbio',
                'mode': 'rw',
            }
        },
        # working_dir="/app",
        command="pytest /root/tests/test_sample_types.py",
        # command='python -c "import unionbio"',
        remove=True,
        stdout=True,
        stderr=True
    )
    
    # Print the output of the pytest command
    print(container)
    # except docker.errors.ContainerError as e:
    #     # Handle errors that occur during container run
    #     print(f"Container error: {e}")
    #     print(e.stderr.decode())
    # except docker.errors.ImageNotFound as e:
    #     # Handle errors if the image is not found
    #     print(f"Image not found: {e}")
    # except docker.errors.APIError as e:
    #     # Handle general API errors
    #     print(f"API error: {e}")

run_pytest_in_docker()

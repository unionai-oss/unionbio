import docker
from config import proj_rt, main_img_test_fqn


def run_pytest_in_docker():
    client = docker.from_env()

    try:
        # Run the Docker container
        print("\nRunning pytest in Docker container...")
        con_name = "unionbio-main-test-container"
        container = client.containers.run(
            name=con_name,
            image=main_img_test_fqn,
            volumes={
                proj_rt.joinpath("tests"): {
                    "bind": "/root/tests",
                    "mode": "rw",
                },
                # Mount unionbio so it's discoverable in PYTHONPATH
                proj_rt.joinpath("src/unionbio"): {
                    "bind": "/root/unionbio",
                    "mode": "rw",
                },
            },
            command="pytest /root/tests/main/test_utils.py::test_intersect_vcfs",
            stdout=True,
            stderr=True,
        )
        # Print the output of the pytest command
        print(container.decode("utf-8"))
    except docker.errors.ContainerError as e:
        # Handle errors that occur during container run
        print(f"Container error: {e}")
        print(e.container.logs().decode("utf-8"))
    except docker.errors.ImageNotFound as e:
        # Handle errors if the image is not found
        print(f"Image not found: {e}")
    except docker.errors.APIError as e:
        # Handle general API errors
        print(f"API error: {e}")
    finally:
        # Ensure the container is cleaned up
        con = client.containers.get(con_name)
        if con is not None:
            con.remove(force=True)
            print(f"Container {con.id} removed")


run_pytest_in_docker()

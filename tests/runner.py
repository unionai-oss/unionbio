import re
import docker
from pathlib import Path
from flytekit.image_spec.image_spec import ImageBuildEngine
from .config import main_img_test_fqn, colabfold_img_test_fqn
from unionbio.images import build_scope
from unionbio.config import project_rt, ws_rt


def run_pytest_in_docker(fqn: str, test_prefix: str, rt: str = None, dr: list = None):
    client = docker.from_env()

    try:
        # Run the Docker container
        print("Running pytest in Docker container...")
        con_name = "unionbio-test-container"
        container = client.containers.run(
            name=con_name,
            image=fqn,
            runtime=rt,
            device_requests=dr,
            volumes={
                project_rt.joinpath("tests"): {
                    "bind": "/root/tests",
                    "mode": "rw",
                },
                # Mount unionbio so it's discoverable in PYTHONPATH
                project_rt.joinpath("src/unionbio"): {
                    "bind": "/root/src/unionbio",
                    "mode": "rw",
                },
            },
            environment={"PYTHONPATH": "/root/src:/root"},
            command=f"pytest /root/tests/{test_prefix}",
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


def test_main():
    run_pytest_in_docker(main_img_test_fqn, "main")


def test_colabfold():
    run_pytest_in_docker(
        colabfold_img_test_fqn,
        "colabfold",
        rt="nvidia",
        dr=[
            docker.types.DeviceRequest(
                count=-1,  # -1 means all available GPUs
                capabilities=[["gpu"]],
            )
        ],
    )


## Images ##


# Updates test config with test image FQNs
def update_img_config(config_path: Path, fqns: dict[str, str]):
    with open(config_path, "r") as f:
        cfg_content = f.read()

    for var, tag in fqns.items():
        cfg_content = re.sub(rf"{var} = .+", f'{var} = "{tag}"', cfg_content)

    with open(config_path, "w") as f:
        f.write(cfg_content)


# Updates workspace YAMLs with test image FQNs
def update_ws_config(config_path: Path, fqn: str):
    with open(config_path, "r") as file:
        yaml_data = file.readlines()

    # Manually parsing instead of using yaml library to preserve structure
    lines_out = []
    for line in yaml_data:
        if line.startswith("container_image:"):
            ll = line.split(": ")
            ll[1] = f"{fqn}\n"
            lines_out.append(": ".join(ll))
        else:
            lines_out.append(line)

    with open(config_path, "w") as file:
        file.writelines(lines_out)


# Entrypoint for building test images in scope
def build_test():
    build_specs = []
    fqns = {}

    # Prepare builds
    for img_str in build_scope:
        spec = eval(img_str)
        spec.tag_format = "{spec_hash}-test"
        spec.source_root = project_rt
        ns = spec.with_packages(["pytest"]).with_apt_packages(["screen", "htop"])
        fqn = ns.image_name()
        fqns[f"{img_str}_test_fqn"] = fqn
        build_specs.append(ns)

        # Update workspace config
        try:
            ws_cfg_path = ws_rt.joinpath(img_str.replace("_img", ""), "ws.yaml")
            update_ws_config(ws_cfg_path, fqn)
        except FileNotFoundError:
            print(f"Workspace config not found for {img_str}")

    for spec in build_specs:
        ImageBuildEngine().build(spec)

    update_img_config(Path("tests/config.py"), fqns)

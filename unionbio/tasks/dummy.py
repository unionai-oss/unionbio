from flytekit import ImageSpec, task

img = ImageSpec(
    name="unionbio-test-2",
    platform="linux/amd64",
    python_version="3.11",
    packages=["flytekit"],
    registry="localhost:30000",
)

@task(container_image=img)
def test_unionbio_install() -> str:
    return "UnionBio package installed successfully."
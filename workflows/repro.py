from flytekit import ImageSpec, task, dynamic

repro_img = ImageSpec(
    name="repro",
    # apt_packages=["curl"],
    packages=["flytekitplugins-envd"],
    registry="ghcr.io/pryce-turner",
    # base_image="alpine:latest"
)

@task(container_image=repro_img)
def foo():
    print("foo")

@dynamic
def bar():
    foo()
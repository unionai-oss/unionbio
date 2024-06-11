from pathlib import Path
from images import ImageFactory, main_img, folding_img, parabricks_img
from flytekit.image_spec.image_spec import ImageBuildEngine

def build_test_images():
    factory = ImageFactory()
    whl = factory.built_wheel()
    # factory.config_path = Path("config.py")
    # fqns = {}

    main_img.with_packages([whl, "pytest"])
    main_img.name = f"{main_img.name}_test"
    ImageBuildEngine().build(main_img)

build_test_images()
    # # Prepare builds
    # for img_str in factory.build_scope:
    #     spec = eval(img_str)
    #     spec.name = f"{spec.name}_test"
    #     spec.with_packages([whl, "pytest"])
    #     fqns[f"{img_str}_fqn"] = spec.image_name()

    # factory.update_img_config(fqns)

    # # Build images
    # for img_str in factory.build_scope:
    #     spec = eval(img_str)
    #     ImageBuildEngine().build(spec)

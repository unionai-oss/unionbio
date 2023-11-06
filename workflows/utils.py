from pathlib import Path
from flytekit.configuration import Config
from flytekit.remote import FlyteRemote, FlyteWorkflowExecution


def get_remote(local=None, config_file=None):
    return FlyteRemote(
        config=Config.auto(
            config_file=(
                None if local
                else config_file if config_file is not None
                else str(Path.home() / ".flyte" / "config-sandbox.yaml")
            )
        ),
        default_project="flytesnacks",
        default_domain="development",
    )
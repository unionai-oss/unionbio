from flytekit import kwtypes, ContainerTask
from flytekit.extras.tasks.shell import ShellTask
from flytekit.experimental import eager
from flytekit.configuration import Config
from flytekit.remote import FlyteRemote


test = ShellTask(
    name="testo",
    debug=True,
    script=
    """
    ls -la
    """,
    inputs=kwtypes(),
    output_locs=[]
)

@eager(
    remote=FlyteRemote(
        config=Config.for_sandbox(),
        default_project="flytesnacks",
        default_domain="development",
    )
)
async def wf():
    out = await _bash.testo()
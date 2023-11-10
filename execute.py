from workflows import alignment, scratch
from workflows.utils import get_remote

remote = get_remote()
execution = remote.execute_local_workflow(
    # alignment.alignment_wf,
    scratch.wf,
    inputs={},
)
print(remote.generate_console_url(execution))

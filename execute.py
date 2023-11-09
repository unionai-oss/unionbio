from workflows import alignment
from workflows.utils import get_remote

remote = get_remote()
execution = remote.execute_local_workflow(
    alignment.alignment_wf,
    inputs={},
)
print(remote.generate_console_url(execution))

from src.workflows.alignment import alignment_wf
from src.tasks.utils import get_remote

remote = get_remote()
execution = remote.execute_local_workflow(
    alignment_wf,
    # scratch.wf,
    inputs={},
)
print(remote.generate_console_url(execution))

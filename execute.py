from src.workflows.simple_alignment import simple_alignment_wf
from src.tasks.utils import get_remote

remote = get_remote()
execution = remote.execute_local_workflow(
    simple_alignment_wf,
    inputs={},
)
print(remote.generate_console_url(execution))

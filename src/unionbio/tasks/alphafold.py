import os
from flytekit import task, workflow, Resources
from flytekit.extras.tasks.shell import subproc_execute
from union.actor import ActorEnvironment
from unionbio.config import alphafold_img_fqn

actor = ActorEnvironment(
    name="af-actor",
    replica_count=1,
    parallelism=1,
    backlog_length=10,
    ttl_seconds=600,
    requests=Resources(
        cpu="8",
        mem="8Gi",
        ephemeral_storage="600Gi",
        gpu="1",
    ),
    container_image=alphafold_img_fqn,
)

@actor.task
def dl_dbs(
    script_loc: str="/app/alphafold/scripts/download_all_data.sh",
    output_loc: str="/mnt/af_dbs",
    reduced: bool=False
    ) -> str:
    os.makedirs(output_loc, exist_ok=True)
    dl_cmd = [
        script_loc,
        output_loc,
    ]
    if reduced:
        dl_cmd.append("reduced_dbs")
    subproc_execute(command=dl_cmd)
    return "\n".join(os.listdir(output_loc))


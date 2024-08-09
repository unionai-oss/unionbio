from unionbio.workflows.simple_variant_calling import calling_wf
from flytekit.remote import FlyteRemote
from flytekit.configuration import Config

remote = FlyteRemote(
    config=Config.for_sandbox(),
    default_project="flytesnacks",
    default_domain="development"
)

# wf = remote.register_script(
#     calling_wf,
#     version="v2",
#     source_path="../",
#     module_name="calling_wf",
# )
# fetched_wf = remote.fetch_workflow(name="workflows.simple_variant_calling.calling_wf")
fet_lp = remote.fetch_launch_plan(name="workflows.simple_variant_calling.calling_wf")
execution = remote.execute(fet_lp,
    inputs={
        "seq_dir": "s3://my-s3-bucket/my-data/sequences",
        "ref_path": "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
    }
)
# execution = remote.execute_local_workflow(
#     # calling_wf,
#     # fetched_wf,
#     inputs={
#         "seq_dir": "s3://my-s3-bucket/my-data/sequences",
#         "ref_path": "s3://my-s3-bucket/my-data/refs/GRCh38_short.fasta"
#     },
# )
print(remote.generate_console_url(execution))

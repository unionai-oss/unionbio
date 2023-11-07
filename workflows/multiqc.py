# from flytekit import kwtypes, task, workflow, ImageSpec, Resources, current_context
# from flytekit.extras.tasks.shell import OutputLocation, ShellTask
# from flytekit.types.file import FlyteFile
# from flytekit.types.directory import FlyteDirectory

# multiqc_image_spec = ImageSpec(
#     name="multiqc",
#     packages=["multiqc"],
#     registry="localhost:30000",
#     base_image='ghcr.io/pryce-turner/variant-discovery:latest'
# )

# multiqc = ShellTask(
#     name="multiqc",
#     debug=True,
#     script=
#     """
#     multiqc {inputs.report_dir} -n {outputs.o}
#     """,
#     inputs=kwtypes(report_dir=FlyteDirectory),
#     output_locs=[OutputLocation(var="o", var_type=FlyteFile, location='/root/multiqc_report.html')],
#     container_image=multiqc_image_spec
# )

# @task(disable_deck=False)
# def render_multiqc(report: FlyteFile):
#     report_html = open(report, 'r').read()
#     current_context().default_deck.append(report_html)

# @workflow
# def multiqc_wf():
#     rep = multiqc(report_dir='s3://my-s3-bucket/data/0n/fb8380196a7e541e5aed-n0-0/3564717cf18a70246792f26354002d68')
#     render_multiqc(report=rep)
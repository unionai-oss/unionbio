from flytekit import TaskMetadata, kwtypes
from flytekit.extras.tasks.shell import OutputLocation, ShellTask
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

from unionbio.config import main_img_fqn

"""
Produce a base quality score recalibration report from a deduped alignment file.

Args:
    rfn (str): The name of the BQSR report file.
    ddal (FlyteFile): An alignment file containing deduped reads.
    ref_fn (str): The basename of the reference fasta.
    ref_dir (FlyteDirectory): The directory containing the reference fasta, index and dict.
    sites (FlyteFile): A VCF file containing known variant sites.

Returns:
    bqsr (FlyteFile): A BQSR report file.
"""
base_recal = ShellTask(
    name="base_recalibrator",
    debug=True,
    metadata=TaskMetadata(retries=3, cache=True, cache_version="1"),
    script="""
    mkdir /tmp/recal
    "java" \
    "-jar" \
    "/usr/local/bin/gatk" \
    "BaseRecalibrator" \
    --input {inputs.ddal} \
    --output {outputs.bqsr} \
    --reference {inputs.ref_dir}/{inputs.ref_fn} \
    --known-sites {inputs.sites} 
    """,
    inputs=kwtypes(
        rfn=str, ddal=FlyteFile, ref_fn=str, ref_dir=FlyteDirectory, sites=FlyteFile
    ),
    output_locs=[
        OutputLocation(
            var="bqsr", var_type=FlyteFile, location="/tmp/recal/{inputs.rfn}"
        )
    ],
    container_image=main_img_fqn,
)

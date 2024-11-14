from flytekit import workflow
from flytekit.types.file import FlyteFile
from unionbio.tasks.colabfold import (
    sync_dbs,
    sync_mmcif,
    cf_search,
    af_predict,
    visualize,
)


@workflow
def cf_wf(db_uris: list[str], mmcif_uri: str, fasta_uri: FlyteFile) -> FlyteFile:
    db_path = sync_dbs(uris=db_uris)
    mmcif_path = sync_mmcif(uri=mmcif_uri)
    hitfile, msa = cf_search(seq=fasta_uri, db_path=db_path)
    af = af_predict(hitfile=hitfile, msa=msa, mmcif_loc=mmcif_path)
    plot = visualize(af_res=af)
    return plot

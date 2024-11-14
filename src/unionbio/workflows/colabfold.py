from flytekit import workflow
from flytekit.types.file import FlyteFile
from unionbio.tasks.colabfold import sync_dbs, sync_mmcif, cf_search, af_predict, visualize

@workflow
def cf_wf() -> FlyteFile:
    db_path = sync_dbs(uris=[
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/cf_envdb/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/pdb100/",
        "gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/uniref30/",
    ])
    mmcif_path = sync_mmcif(uri="gs://opta-gcp-dogfood-gcp/bio-assets/colabfold/mmcif_tar/")
    hitfile, msa = cf_search(seq="gs://opta-gcp-dogfood-gcp/bio-assets/fastas/rcsb_pdb_2HHB.fasta", db_path=db_path)
    af = af_predict(hitfile=hitfile, msa=msa, mmcif_loc=mmcif_path)
    plot = visualize(af_res=af)
    return plot
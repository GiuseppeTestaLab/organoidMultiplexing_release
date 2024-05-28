import scanpy as sc


def PurgeAdata(adata):
	_adata = adata.copy()
	for i in [i for i in list(_adata.uns.keys()) if "_colors" not in i]:
		del _adata.uns[i]
	if "leiden" in _adata.obs.columns:
		del _adata.obs["leiden"]
	if "louvain" in _adata.obs.columns:
		del _adata.obs["louvain"]
	if "highly_variable" in _adata.var.columns:
		del _adata.var["highly_variable"]
	del _adata.varm
	del _adata.obsp 
	del _adata.obsm
	return _adata

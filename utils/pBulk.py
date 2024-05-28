import scanpy as sc
import pandas as pd
import random
import itertools
import os
from scipy.sparse import  csr_matrix, issparse


def pBulk(adata, PseudooReplicates_per_group , groupingCovariate, normalizepre=False,normalizepost=False,layer=None, seed=1):
	_adata = adata.copy()
	# Select layer
	if layer is None:
		print("No layer selected, counts in .X will be aggregated")
		_adata.layers["aggrLayer"] = _adata.X.copy()
	elif layer is not None:
		print("Counts in {} will be aggregated".format(layer))
		_adata.layers["aggrLayer"] = _adata.layers[layer].copy()
	if not issparse(_adata.layers["aggrLayer"]):
		_adata.layers["aggrLayer"] = csr_matrix(_adata.layers["aggrLayer"])
	if normalizepre:
		sc.pp.normalize_total(_adata, target_sum=20000, layer="aggrLayer")
	print("Pbulking with target of "+str(PseudooReplicates_per_group)+" PRs will result in following cells per PR")
	print(_adata.obs.groupby(groupingCovariate).size() / PseudooReplicates_per_group)
	
	total = pd.DataFrame(index = _adata.var_names)
	total_metadata = pd.DataFrame(index= ["@".join(Rep) for Rep in  list(itertools.product(_adata.obs[groupingCovariate].cat.categories.tolist(), [str(r)  for r in list(range(PseudooReplicates_per_group))]))])
	
	metaCellsDict = {}
	for group in _adata.obs[groupingCovariate].cat.categories:
		groupAdata = _adata[_adata.obs[groupingCovariate] == group]
		
		group_cells = groupAdata.obs_names.tolist()
		random.Random(seed).shuffle(group_cells)
		
		# if os.path.isfile(homeDir+"/data/resources/{}.selectedBCs.txt".format(nb_fname)):
		# 	metaCellslist = pd.read_csv(homeDir+"/data/resources/{}.selectedBCs.txt".format(nb_fname), sep="\t", index_col=0)
		# 	metaCellslist = metaCellslist.loc[metaCellslist.index.str.startswith(str(group)+"@"),"Barcode"].apply(lambda row: row.split(",")).values.tolist()
		# else:
		# 	metaCellslist=[group_cells[i::PseudooReplicates_per_group] for i in range(PseudooReplicates_per_group)]
		metaCellslist=[group_cells[i::PseudooReplicates_per_group] for i in range(PseudooReplicates_per_group)]
		for partition in enumerate(metaCellslist):
			
			total['{}@{}'.format(group, partition[0])] = _adata[partition[1]].layers["aggrLayer"].sum(axis = 0).A1
			metaCellsDict['{}@{}'.format(group, partition[0])] = partition[1]
		
			total_metadata.loc['{}@{}'.format(group, partition[0]),groupingCovariate] = group
			total_metadata.loc['{}@{}'.format(group, partition[0]),"pseudoreplicate"] = partition[0]
			total_metadata.loc['{}@{}'.format(group, partition[0]),"number_of_cell"] = int(len(partition[1]))
			total_metadata.loc['{}@{}'.format(group, partition[0]),"midexBCs"] = ",".join(partition[1])
		#With this line we propagate - whenever possible -  obs to aggregated pReplicate
		for obsMD in [obsMD for obsMD in groupAdata.obs.columns.tolist() if len(groupAdata.obs[obsMD].unique()) == 1 and obsMD != groupingCovariate]:
			total_metadata.loc[["@".join(l) for l in list(itertools.product([group],[str(r)  for r in list(range(PseudooReplicates_per_group))]))], obsMD ] = _adata.obs.loc[_adata.obs[groupingCovariate] == group,obsMD][0]
			
			
	total_metadata = total_metadata.dropna(axis = 1)
	
	
	
	adata_bulk = sc.AnnData(total.transpose()).copy()
	adata_bulk.var = _adata.var.copy()
	adata_bulk.obs = pd.concat([adata_bulk.obs, total_metadata], axis = 1)

	adata_bulk.obs[groupingCovariate] =adata_bulk.obs[groupingCovariate].astype("category", ord)
	try:
		cmap = dict(zip(adata.obs[groupingCovariate].cat.categories.tolist(), list(adata.uns["{}_colors".format(groupingCovariate)])))
		cmap
		adata_bulk.obs[groupingCovariate] = adata_bulk.obs[groupingCovariate].cat.reorder_categories(list(cmap.keys()))
		adata_bulk.uns["{}_colors".format(groupingCovariate)] = list(cmap.values())
	except:
		print("No color mapping was possible for {}".format(groupingCovariate))
	
	adata_bulk.uns["pSamples"] = metaCellsDict
	return adata_bulk



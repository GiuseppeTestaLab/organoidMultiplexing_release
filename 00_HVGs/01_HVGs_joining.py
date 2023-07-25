#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import anndata as ad
import yaml
import socket
import os
import warnings


# In[2]:



outdir = "../data/output"
if not os.path.exists(outdir):
   # Create a new directory because it does not exist
   os.makedirs(outdir)

with open("../data/resources/iPSC_lines_map.yaml", 'r') as f:
    iPSC_lines_map = yaml.load(f, Loader=yaml.FullLoader)["lines"]

indir = "../data"


# In[3]:


UpD50Path = "../data/Sample_S20272_157/filtered_feature_bc_matrix/"
DownD50Path = "../data/Sample_S20273_158/filtered_feature_bc_matrix/"


UpD100_1Path = "../data/Sample_S20812_258/filtered_feature_bc_matrix/"
UpD100_2Path =  "../data/Sample_S20813_259/filtered_feature_bc_matrix/"
DownD100Path = "../data/Sample_S31807_MET6_GEX/filtered_feature_bc_matrix/"


UpD300Path = "../data/Sample_S33846_C_GEX/filtered_feature_bc_matrix/"
DownD250Path = "../data/Sample_S20814_260/filtered_feature_bc_matrix/"


# ### D50

# In[4]:


#D50 HVG <------------------------------------------------------------------------------
#D50 HVG <------------------------------------------------------------------------------


#UpD50-----------------------------------------
#UpD50-----------------------------------------

DSname="UpD50"
DSnameDirName="Sample_S20272_157"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 9, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 6.5, :]

adata = adata[adata.obs.log1p_total_counts_mt < 6.5, :]
adata = adata[adata.obs.log1p_total_counts_mt > 1.5, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 9, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 4.5, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
UpD50 = adata.copy()


#DownD50-----------------------------------------
#DownD50-----------------------------------------

DSname="DownD50"
DSnameDirName="Sample_S20273_158"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 9, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 6.5, :]

adata = adata[adata.obs.log1p_total_counts_mt < 6.5, :]
adata = adata[adata.obs.log1p_total_counts_mt > 1.5, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 9, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 4.5, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
DownD50 = adata.copy()


del adata


UpD50.obs['dataset'] = "UpD50"
UpD50.obs_names = [i + "_" + j for i, j in zip(UpD50.obs_names.tolist(), UpD50.obs["dataset"].tolist())]
DownD50.obs['dataset'] = "DownD50"
DownD50.obs_names = [i + "_" + j for i, j in zip(DownD50.obs_names.tolist(), DownD50.obs["dataset"].tolist())]
#Merge datasets
adata = ad.concat([UpD50,DownD50], merge="same")


sc.pp.normalize_total(adata, target_sum=2e4, layers = "all")
sc.pp.log1p(adata)
#scv.pp.filter_genes(adata, min_shared_counts = 30)
#HVG detection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0.5, batch_key = "dataset")
#adata = adata[:, adata.var.highly_variable]
HVGd50 = adata.var[adata.var.highly_variable_nbatches >= 2].index.tolist()
HVGd50Union = adata.var[adata.var.highly_variable_nbatches >= 1].index.tolist()


# ### D100

# In[5]:


#D100 HVG <--------------------------
#D100 HVG <--------------------------
#D100 HVG <--------------------------

#UpD100_1-----------------------------------------
#UpD100_1-----------------------------------------

DSname="UpD100_1"
DSnameDirName="Sample_S20812_258"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 9, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 6.5, :]

adata = adata[adata.obs.log1p_total_counts_mt < 6.5, :]
adata = adata[adata.obs.log1p_total_counts_mt > 1.5, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 9, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 4.5, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
UpD100_1 = adata.copy()


#UpD100_2-----------------------------------------
#UpD100_2-----------------------------------------

DSname="UpD100_2"
DSnameDirName="Sample_S20813_259"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 9, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 6.5, :]

adata = adata[adata.obs.log1p_total_counts_mt < 6.5, :]
adata = adata[adata.obs.log1p_total_counts_mt > 1.5, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 9, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 4.5, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
UpD100_2 = adata.copy()


#DownD100-----------------------------------------
#DownD100-----------------------------------------

DSname="DownD100"
DSnameDirName="Sample_S31807_MET6"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 8, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 6, :]

adata = adata[adata.obs.log1p_total_counts_mt < 7, :]
adata = adata[adata.obs.log1p_total_counts_mt > 5.5, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 5.5, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 2, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
DownD100 = adata.copy()


del adata


UpD100_1.obs['dataset'] = "UpD100_1"
UpD100_1.obs_names = [i + "_" + j for i, j in zip(UpD100_1.obs_names.tolist(), UpD100_1.obs["dataset"].tolist())]

UpD100_2.obs['dataset'] = "UpD100_2"
UpD100_2.obs_names = [i + "_" + j for i, j in zip(UpD100_2.obs_names.tolist(), UpD100_2.obs["dataset"].tolist())]

DownD100.obs['dataset'] = "DownD100"
DownD100.obs_names = [i + "_" + j for i, j in zip(DownD100.obs_names.tolist(), DownD100.obs["dataset"].tolist())]
#Merge datasets
adata = ad.concat([UpD100_1,UpD100_2, DownD100], merge="same")



sc.pp.normalize_total(adata, target_sum=2e4, layers = "all")
sc.pp.log1p(adata)
#scv.pp.filter_genes(adata, min_shared_counts = 30)
#HVG detection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0.5, batch_key = "dataset")
#adata = adata[:, adata.var.highly_variable]
HVGd100 = adata.var[adata.var.highly_variable_nbatches >= 2].index.tolist()
HVGd100Union = adata.var[adata.var.highly_variable_nbatches >= 1].index.tolist()


# ### D250

# In[6]:


#D250 HVG <--------------------------
#D250 HVG <--outdir----------------------
#D250 HVG <--------------------------


#UpD300-----------------------------------------
#UpD300-----------------------------------------

DSname="UpD300"
DSnameDirName="Sample_S33846_C_GEX"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 8.5, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 5.5, :]

adata = adata[adata.obs.log1p_total_counts_mt < 7.5, :]
adata = adata[adata.obs.log1p_total_counts_mt > 3, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 8, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 4.5, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
UpD300 = adata.copy()


#DownD250-----------------------------------------
#DownD250-----------------------------------------

DSname="DownD250"
DSnameDirName="Sample_S20814_260"
adata=sc.read_10x_mtx(indir+'/'+DSnameDirName+"/filtered_feature_bc_matrix/",var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
cellID = pd.read_csv(indir+'/'+DSnameDirName+'/aggregatedCall/aggregatedCall.tsv', sep ="\t",index_col = 0)
adata.obs['cellID'] = cellID.loc[adata.obs_names, 'Consensus']
adata = adata[adata.obs['cellID'].isin(["3391B","KOLF","MIFF1","809","LowQuality","doublet"])]
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

adata.var['ribo'] = adata.var_names.str.startswith('RP')  # annotate the group of ribosomal genes as 'ribo'
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=True, inplace=True)


adata = adata[adata.obs.log1p_n_genes_by_counts < 9, :]
adata = adata[adata.obs.log1p_n_genes_by_counts > 6.5, :]

adata = adata[adata.obs.log1p_total_counts_mt < 6.5, :]
adata = adata[adata.obs.log1p_total_counts_mt > 1.5, :]

adata = adata[adata.obs.log1p_total_counts_ribo < 9, :]
adata = adata[adata.obs.log1p_total_counts_ribo > 4.5, :]
adata = adata[-adata.obs['cellID'].isin(["LowQuality","doublet"])]
pd.DataFrame(adata.obs_names, columns=["BCs"]).to_csv(outdir+'/'+DSname+'_filteredCells.tsv', sep="\t")
DownD250 = adata.copy()


del adata

UpD300.obs['dataset'] = "UpD300"
UpD300.obs_names = [i + "_" + j for i, j in zip(UpD300.obs_names.tolist(), UpD300.obs["dataset"].tolist())]
DownD250.obs['dataset'] = "DownD250"
DownD250.obs_names = [i + "_" + j for i, j in zip(DownD250.obs_names.tolist(), DownD250.obs["dataset"].tolist())]
#Merge datasets
adata = ad.concat([UpD300,DownD250], merge="same")


sc.pp.normalize_total(adata, target_sum=2e4, layers = "all")
sc.pp.log1p(adata)
#scv.pp.filter_genes(adata, min_shared_counts = 30)
#HVG detection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0.5, batch_key = "dataset")
#adata = adata[:, adata.var.highly_variable]
HVGd250 = adata.var[adata.var.highly_variable_nbatches >= 2].index.tolist()
HVGd250Union = adata.var[adata.var.highly_variable_nbatches >= 1].index.tolist()


HVGoverlap = list(set(HVGd100+HVGd50+HVGd250))
HVUnionsGoverlap = list(set(HVGd250Union+HVGd50Union+HVGd100Union))


# In[7]:




External = list(set(["MKI67","CDC20","PAX6","NES","HOPX",
"DCX","STMN2","NEUROD6","NEUROD2","MAP2","SYN1",
"SLC17A6", "SATB2","BCL11B","TBR1", "GAD2","DLX2","DLX5","SLC32A1","S100B", 
"GFAP", "AQP4", "OLIG1","DCN", "BGN", "MYH3", "GAPDH",
"OLIG1","OLIG2", "OLIG3","PDGFRA",
"HOPX","S100B", "GFAP", "AQP4", "TNC",
"GAD1","GAD2","DLX1","DLX2","SLC32A1","DLX6-AS1",
"CA8","FOXP1",
"SLC17A6","TBR1","SLA",
"NEUROD2","NEUROD6","RELN","NEUROD4",
"GAP43","DCX","STMN2","MAP2","SYN1","NEUROD2","NEUROD6","CUX1",
"SOX2","PAX6","NES","VIM","HES1",
"MKI67","CDC20",]))


HVGoverlap_curated = list(set(HVGoverlap+External))
pd.DataFrame(HVGoverlap_curated, columns = ["HVG"]).to_csv(outdir+"/HVG_list_intersection_Curated.txt", sep="\t")


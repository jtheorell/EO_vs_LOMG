#Here, we are dependent on: 
#https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/cite_scrna_integration_w_totalVI.html
#This relates to the totalVI publication: https://pubmed.ncbi.nlm.nih.gov/33589839/

python3

import tempfile
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import scanpy as sc
import scvi
import torch
from scipy.stats import pearsonr

os.chdir("/Users/jakob.theorell/Labbet/2024/240829_merged_MG_analyses/External/Stockholm/Data")

scvi.settings.seed = 0
sc.set_figure_params(figsize=(4, 4), frameon=False)
torch.set_float32_matmul_precision("high")

allGeneData = pd.read_csv('Anndata_input/geneData.csv', index_col=0)
allCellData = pd.read_csv('Anndata_input/cellData.csv', low_memory = False, index_col=0)
allAdtData = pd.read_csv("Anndata_input/adt_counts_common.csv", index_col=0)

adata_myakot = sc.read_mtx("Anndata_input/gex_counts.mtx")
adata_myakot.var = allGeneData
adata_myakot.obs = allCellData
adata_myakot.obsm["protein_expression"] = allAdtData.values

sc.pp.highly_variable_genes(
    adata_myakot, batch_key="batch", flavor="seurat_v3", n_top_genes=4000
)
adata_myakot.write("1_adata_common_complete.h5ad")
#adata_myakot = sc.read_h5ad('1_adata_common_complete.h5ad')
scvi.model.TOTALVI.setup_anndata(
    adata_myakot, batch_key="batch", protein_expression_obsm_key="protein_expression"
)

model = scvi.model.TOTALVI(adata_myakot, latent_distribution="normal", n_layers_decoder=2)

model.train(use_gpu = False, max_epochs = 200)

x_totalVi = pd.DataFrame(model.get_latent_representation())
os.makedirs("Normalised_protein_models_common_proteins")
x_totalVi.to_csv("Normalised_protein_models_common_proteins/Latent_representation_common_proteins.csv")

#To avoid braking the computer, I have divided the below into 25 runs. 
batch_names = ['Ctrl_1', 'Ctrl_2', 'Pat_1', 'Pat_2']

for x in range(25):
    print(x)
    scvi.settings.seed = x
    total_model = model.get_normalized_expression(transform_batch=batch_names, n_samples=1)
    total_model_prot = pd.DataFrame(total_model[1])
    file_name = "".join(["Normalised_protein_models_common_proteins/Model_", str(x), ".csv"])
    total_model_prot.to_csv(file_name)

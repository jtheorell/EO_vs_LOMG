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

#Here, we add the local directory where the files are stored, in my case
os.chdir("/Users/jakob.theorell/Labbet/2024/240829_merged_MG_analyses/External/Stockholm/Data")

scvi.settings.seed = 0
sc.set_figure_params(figsize=(4, 4), frameon=False)
torch.set_float32_matmul_precision("high")

allGeneData = pd.read_csv('Anndata_input/geneData.csv', index_col=0)
allCellData = pd.read_csv('Anndata_input/cellData.csv', low_memory = False, index_col=0)
allAdtData = pd.read_csv("Anndata_input/adt_counts.csv", index_col=0)

adata_myakot = sc.read_mtx("Anndata_input/gex_counts.mtx")
adata_myakot.var = allGeneData
adata_myakot.obs = allCellData
adata_myakot.obsm["protein_expression"] = allAdtData.values

sc.pp.highly_variable_genes(
    adata_myakot, batch_key="batch", flavor="seurat_v3", n_top_genes=4000
)
adata_myakot.write("1_adata_complete.h5ad")
#adata_myakot = sc.read_h5ad('1_adata_complete.h5ad')

scvi.model.TOTALVI.setup_anndata(
    adata_myakot, batch_key="batch", protein_expression_obsm_key="protein_expression"
)

model = scvi.model.TOTALVI(adata_myakot, latent_distribution="normal", n_layers_decoder=2)

model.train(use_gpu = False, max_epochs = 200)

#These are saved for now, code-wise, but the data is not, as it is the average of the normalised models that we will use. 
x_totalVi = pd.DataFrame(model.get_latent_representation())
os.makedirs("Normalised_protein_models")
x_totalVi.to_csv("Normalised_protein_models/Latent_representation.csv")

#We batch this to handle it without breaking the computer. 
batch_names = ['Ctrl_1', 'Ctrl_2', 'Pat_1', 'Pat_2']

for x in range(25):
    print(x)
    scvi.settings.seed = x
    total_model = model.get_normalized_expression(transform_batch=batch_names, n_samples=1)
    total_model_prot = pd.DataFrame(total_model[1])
    file_name = "".join(["Normalised_protein_models/Model_", str(x), ".csv"])
    total_model_prot.to_csv(file_name)
#This worked, even if for some reason, the 24:th iteration did not. I had however
#already saved the 25th, so that meant a total of 24 runs that we can average. 
#And after this, we go back to R. 

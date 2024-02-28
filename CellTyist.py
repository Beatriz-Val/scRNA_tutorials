"""
Created on 15/02/2024, 10:01, 2024

@autor: Beatriz Val 68818 MCBBI

"""

# Import libraries
import pandas as pd
import scanpy as sc  # single-cell package
import celltypist  # for the classification
import anndata2ri  # conversion
from rpy2.robjects import r  # running R in python
anndata2ri.activate()  # activate automatic SCE<->AnnData conversion

# Open seurat object in R
r('obj = readRDS("merged_integrated_annotated.rds")')


# convert Seurat to SCE
adata_gn_l = r('Seurat::as.SingleCellExperiment(obj)')


# Download models
celltypist.models.download_models(model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])


sc.pp.normalize_total(adata_gn_l, target_sum=10000)
sc.pp.log1p(adata_gn_l)

pred_dic = celltypist.annotate(adata_gn_l, model='Immune_All_Low.pkl', majority_voting=True)

print(pred_dic)

pred_dic.predicted_labels.to_csv('pred_dic.csv')

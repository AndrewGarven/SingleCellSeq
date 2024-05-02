import scanpy as sc
from harmony import harmonize


adata = sc.read_h5ad('/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/Harmony/Lai/Input/input_harmony.h5ad')
Z = harmonize(adata.obsm['X_pca'], adata.obs, batch_key = 'sample')
adata.obsm['X_harmony'] = Z

sc.write("/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/Harmony/Lai/output/output_harmony.h5ad", adata)

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
seurat_obj<- YDL
# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
seurat_obj$tSNE_1 <- seurat_obj@reductions$tsne@cell.embeddings[,1]
seurat_obj$tSNE_2 <- seurat_obj@reductions$tsne@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')
# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)





#conda activate velocity

conda activate velocyto
python
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import matplotlib as mpl

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
mpl.use('Agg')


# load sparse matrix:
X = io.mmread("counts.mtx")

# create anndata object
adata = anndata.AnnData(
  X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")

# load gene names:
with open("gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
adata.obsm['X_tsne'] = np.vstack((adata.obs['tSNE_1'].to_numpy(), adata.obs['tSNE_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['seurat_clusters'], frameon=False, save="True.pdf")

# save dataset as anndata format
adata.write('my_data.h5ad')

# reload dataset
adata = sc.read_h5ad('my_data.h5ad')


ldata1=scv.read_loom("/data3/yudonglin/ren/HRR016237/velocyto/HRR016237.loom",validate=False)
ldata2=scv.read_loom("/data3/yudonglin/ren/HRR016238/velocyto/HRR016238.loom",validate=False)
ldata3=scv.read_loom("/data3/yudonglin/ren/HRR016239/velocyto/HRR016239.loom",validate=False)
ldata4=scv.read_loom("/data3/yudonglin/ren/HRR016240/velocyto/HRR016240.loom",validate=False)



# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = ["liver metastasized tumor from rectal cancer_"+bc[0:len(bc)-1] for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = ["primary rectal cancer_"+bc[0:len(bc)-1] for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = ["adjacent non-tumor tissue of hepatocellular carcinoma_"+bc[0:len(bc)-1] for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = ["liver tumor_"+bc[0:len(bc)-1] for bc in barcodes]
ldata4.obs.index = barcodes


# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()



# concatenate the three loom
ldata = ldata1.concatenate([ldata1, ldata2, ldata3, ldata4])         
# concatenate the three loom
#ldata = ldata1.concatenate([ldata1, ldata2, ldata3, ldata4, ldata5])
# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
# plot umap to check
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='', save='_celltypes.pdf')

sc.pl.tsne(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='', save='tsne.pdf')

scv.pl.proportions(adata, groupby='seurat_clusters',save="proportion.pdf")


# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)


scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding_umap.pdf')
scv.pl.velocity_embedding(adata, basis='tsne', frameon=False, save='embedding_tsne.pdf')


scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save='embedding_grid_umap.pdf', title='', scale=0.25)

scv.pl.velocity_embedding_grid(adata, basis='tsne', color='seurat_clusters', save='embedding_grid_tsne.pdf', title='', scale=0.25)


scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', save='embedding_stream_umap.pdf', title='')
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='seurat_clusters', save='embedding_stream_tsne.pdf', title='')

#可以调节长宽比，不过最后只能输出图片，而不能是PDF
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', figsize =(10, 10),save='embedding_stream_umap1.pdf', title='')
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='seurat_clusters', figsize =(10, 10),save='embedding_stream_tsne1.pdf', title='')



scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltype', save='embedding_grid_umap_celltype.pdf', title='', scale=0.25)

scv.pl.velocity_embedding_grid(adata, basis='tsne', color='celltype', save='embedding_grid_tsne_celltype.pdf', title='', scale=0.25)


scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype', save='embedding_stream_umap_celltype.pdf', title='')
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='celltype', save='embedding_stream_tsne_celltype.pdf', title='')

#可以调节长宽比，不过最后只能输出图片，而不能是PDF
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype', figsize =(10, 10),save='embedding_stream_umap1——celltype.pdf', title='')
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='celltype', figsize =(10, 10),save='embedding_stream_tsne1——celltype.pdf', title='')




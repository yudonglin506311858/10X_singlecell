
conda create -n velocity 
conda activate velocity 
conda install R
#然后安装需要的库：
conda install -c conda-forge openmp
conda install -c anaconda boost
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
conda install r-devtools
conda install r-pkgdown
conda install r-ragg


install.packages("devtools")
library(devtools)
install_github("velocyto-team/velocyto.R")
library(velocyto.R)







conda activate velocity
export OPENBLAS_NUM_THREADS=1

# #导入包
# set.seed(123)
# library(Seurat)
# library(velocyto.R)
# library(pagoda2)
# #导入包
# library(Seurat)
# library(tidyverse)
# library(patchwork)
# library(SCENIC)
# library(plyr)
# library(permute)
# library(data.table)
# library(SCopeLoomR)
# library(velocyto.R)
# library(pagoda2)
# library(devtools)
# install_github("velocyto-team/velocyto.R")
# BiocManager::install("pcaMethods")
# BiocManager::install("collapse")
# BiocManager::install("collapse")
# library(devtools)
# install_github("velocyto-team/velocyto.R")
# #导入包
# set.seed(123)
# library(Seurat)
# library(velocyto.R)
# library('hdf5r')
# library(pcaMethods)
# library(pcaMethods)
# install_github("velocyto-team/velocyto.R")
# #导入包
# set.seed(123)
# library(Seurat)
# library(velocyto.R)
# library(pagoda2)
# #导入包
# library(Seurat)
# library(tidyverse)
# library(patchwork)
# library(SCENIC)
# library(plyr)
# library(permute)
# library(data.table)
# library(pagoda2)
# library(SCopeLoomR)
# library(velocyto.R)
# library(devtools)
# install_github("velocyto-team/velocyto.R")
# #导入包
# set.seed(123)
# library(Seurat)
# library(velocyto.R)
# devtools::install_github('kharchenkolab/pagoda2')
# Installing Linux dependencies
#导入包
set.seed(123)
library(Seurat)
library(velocyto.R)
set.seed(123)
library(Seurat)
library(velocyto.R)
library(pagoda2)
#导入包
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(velocyto.R)
library(pagoda2)
set.seed(123)
library(Seurat)
library(velocyto.R)
library(pagoda2)
#导入包
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(velocyto.R)
library(pagoda2)
set.seed(123)
library(Seurat)
library(velocyto.R)
library(pagoda2)
#导入包
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(velocyto.R)
library(pagoda2)
set.seed(123)
library(Seurat)
library(velocyto.R)
library(pagoda2)
#导入包
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(velocyto.R)
library(pagoda2)
setwd("/data3/yudonglin/ren")
#YDL<- readRDS("/data3/yudonglin/ren/cell_matastasis.RDS")

#liver metastasized tumor from rectal cancer
x1 <-velocyto.R::read.loom.matrices("/data3/yudonglin/ren/HRR016237/velocyto/HRR016237.loom")
splice<-x1$spliced
unsplice<-x1$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste("liver metastasized tumor from rectal cancer_",substring(colnames(splice),11,26),"-1",sep="")
colnames(unsplice) <- paste("liver metastasized tumor from rectal cancer_",substring(colnames(unsplice),11,26),"-1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x1$spliced<-splice
x1$unspliced<-unsplice
table(YDL@meta.data$orig.ident)
#primary rectal cancer
x2 <-velocyto.R::read.loom.matrices(file = "/data3/yudonglin/ren/HRR016238/velocyto/HRR016238.loom", engine = "hdf5r")
splice<-x2$spliced
unsplice<-x2$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste("primary rectal cancer_",substring(colnames(splice),11,26),"-1",sep="")
colnames(unsplice) <- paste("primary rectal cancer_",substring(colnames(unsplice),11,26),"-1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x2$spliced<-splice
x2$unspliced<-unsplice
table(YDL@meta.data$orig.ident)
#adjacent non-tumor tissue of hepatocellular carcinoma
x3 <-velocyto.R::read.loom.matrices(file = "/data3/yudonglin/ren/HRR016239/velocyto/HRR016239.loom", engine = "hdf5r")
splice<-x3$spliced
unsplice<-x3$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste("adjacent non-tumor tissue of hepatocellular carcinoma_",substring(colnames(splice),11,26),"-1",sep="")
colnames(unsplice) <- paste("adjacent non-tumor tissue of hepatocellular carcinoma_",substring(colnames(unsplice),11,26),"-1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x3$spliced<-splice
x3$unspliced<-unsplice
table(YDL@meta.data$orig.ident)
#liver tumor
x4 <-velocyto.R::read.loom.matrices(file = "/data3/yudonglin/ren/HRR016240/velocyto/HRR016240.loom", engine = "hdf5r")
splice<-x4$spliced
unsplice<-x4$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste("liver tumor_",substring(colnames(splice),11,26),"-1",sep="")
colnames(unsplice) <- paste("liver tumor_",substring(colnames(unsplice),11,26),"-1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x4$spliced<-splice
x4$unspliced<-unsplice
spliced <- cbind(x1[["spliced"]], x2[["spliced"]],x3[["spliced"]], x4[["spliced"]])
unspliced <- cbind(x1[["unspliced"]], x2[["unspliced"]],x3[["unspliced"]], x4[["unspliced"]])
emat <- spliced
nmat <- unspliced
seurat.object<-YDL
emb <- seurat.object@reductions$tsne@cell.embeddings
# Estimate the cell-cell distances
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$tsne@cell.embeddings)))
fit.quantile <- 0.02
# Main velocity estimation
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=48)
# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
library("Seurat")
library("ggplot2")
gg <- TSNEPlot(seurat.object)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)
p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=60,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                     do.par=F,cell.border.alpha = 0.1,
                                     n.cores=48,main="RNA Velocity")
p2<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                                   cell.colors = ac(colors, alpha = 0.5),
                                   cex = 0.8, arrow.scale = 25, show.grid.flow = T,
                                   min.grid.cell.mass = 0.5, grid.n = 40,
                                   arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5)
p3<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                                   cell.colors = ac(colors, alpha = 0.5),
                                   cex = 0.8, arrow.scale = 50, show.grid.flow = T,
                                   min.grid.cell.mass = 0.5, grid.n = 40,
                                   arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5)
dev.off()
savehistory("code.txt")


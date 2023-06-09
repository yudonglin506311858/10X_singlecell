

YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
##存储marker

write.csv(top10 ,file="allmarker_cluster.csv")
write.csv(YDL.markers,file="allmarker_cluster.csv")
write.csv(YDL.markers,file="allmarker_celltype.csv")
YDL.markers<-read.csv("allmarker_YDL.csv")
#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_YDL.pdf",width=12,height=9)
DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()



library(monocle)
#准备monocle分析需要的文件
monocle.matrix=as.matrix(YDL@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(YDL@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(YDL.markers,file="monocleMarkers.txt",sep="\t",row.names=F,quote=F)


#设置工作目录
monocle.matrix=read.table("monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("monocleMarkers.txt",sep="\t",header=T,check.names=F)

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表表
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)


#给其中一列数据重命名
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
# saveRDS(cds,"cds.rds")
# rm(list = ls())  
# marker<-read.csv("allmarker.csv")
# cds<-readRDS("cds.rds")
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds)
pdf(file="cluster.trajectory_SINGLET.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~Cluster,nrow=3,ncol = 3)
plot_cell_trajectory(cds,color_by="orig.ident")+facet_wrap(~orig.ident,nrow=2,ncol = 5)
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~orig.ident,nrow=2,ncol = 3)

plot_cell_trajectory(cds, color_by = "YDL@active.ident")
plot_cell_trajectory(cds, color_by="Pseudotime", show_backbone=FALSE)
# 可以很明显看到细胞的发育轨迹 
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1)
dev.off()

saveRDS(cds,"cds_SINGLET.rds")





#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

YDL<-readRDS("cell_matastasis.RDS")
DimPlot(YDL, reduction = "umap",pt.size = 1.5,label = T)
###Pseudotime monocle3
cds <- as.cell_data_set(YDL)
cds <- cluster_cells(cds)
#head(pData(cds))
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


stem1<-rownames(pData(cds)[which(pData(cds)$ident %in% c('B cells')),])
cds <- order_cells(cds, root_cells = stem1)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

stem2<-rownames(pData(cds)[which(pData(cds)$ident %in% c('fibroblasts')),])
cds <- order_cells(cds, root_cells = stem2)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)



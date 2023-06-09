Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
setwd("/data3/yudonglin/ren/10X")
#批量读取10X单细胞转录组数据文件夹：
folders=list.files('./',pattern='[ar]$')
folders
library(Seurat)
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
YDL <- merge(sceList[[1]], 
             y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
             add.cell.ids = folders, 
             project = "rectal cancer")
YDL
head(YDL)
table(YDL$orig.ident)
rm(sceList)
#https://zhuanlan.zhihu.com/p/121439929

set.seed(123)
getwd()
##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中，mit-开头的为线粒体基因
YDL[["percent.mt"]] <- PercentageFeatureSet(object = YDL, pattern = "^MT-")
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)


head(YDL@meta.data)
table(YDL@meta.data$orig.ident)
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)


##过滤细胞：根据上面小提琴图中基因数"nFeature_RNA"和线粒体数"percent.mt"，分别设置过滤参数，这里基因数 200-4000，线粒体百分比为小于 5%，保留gene数大于200小于2500的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于5%的细胞,为了过滤掉死细胞等低质量的细胞数据。
YDL <- subset(x = YDL, subset = nFeature_RNA > 200 &nFeature_RNA < 7000 &nCount_RNA<30000 &nCount_RNA>1000 & percent.mt < 15)    #对数据进行过滤
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Idents(YDL))
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)


head(YDL@meta.data)

YDL<-readRDS("cell_matastasis.RDS")
saveRDS(YDL,"cell_matastasis.RDS")

#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 50,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL,ndims = 50)#选择top20个PC
pcSelect=40
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution = 0.3)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包，这里不做介绍，两种方法都可以使用，但不要混用，如果混用，后面的结算结果会将先前的聚类覆盖掉，只能保留一个。

##这里采用基于TSNE的聚类方法。
YDL <- RunTSNE(object = YDL, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_proe.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化
#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")
dev.off()
write.table(YDL$seurat_clusters,file="tsneCluster_proe.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
YDL <- RunUMAP(object = YDL, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__proe.pdf",width=6.5,height=6)
UMAPPlot(object = YDL, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(object=YDL,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

YDL<- CellCycleScoring(object = YDL, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = YDL@meta.data)
pdf(file="CellCycle_proe.pdf",width=6.5,height=6)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()



YDL[["celltype"]] <- Idents(YDL)
table(YDL@meta.data$celltype)#

#细胞比例变化
library("Seurat")
library("ggplot2")
gg <- TSNEPlot(YDL)
col<-ggplot_build(gg)$data
col<-as.data.frame(col)
table(col$colour)
table(YDL$seurat_clusters)
a<-as.data.frame(table(col$colour))
b<-as.data.frame(table(YDL$orig.ident))
c<-merge(a,b,by = "Freq")
c<-c[order(c$Freq,decreasing = T),]
# 计算细胞比例
Idents(YDL)<-YDL@meta.data$seurat_clusters
prop.table(table(Idents(YDL)))
#pdf("细胞比例——柱状图.pdf")
#col <- c("#F8766D","#CD9600","#7CAE00","#8CAB00","#24B700","#00BE70","#00C1AB","#00BBDA","#00ACFC","#F8766D","#D575FE","#F962DD","#FF65AC")
col<- as.character(c$Var1.x)
## 绘制堆叠条形图
library("ggplot2")
cell.prop<-as.data.frame(prop.table(table(Idents(YDL), YDL$orig.ident)))
colnames(cell.prop)<-c("cluster","sample","proportion")
ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+ theme(axis.text.x = element_text(angle=45, hjust=1)) 

#https://www.jianshu.com/p/56c1dcf3d9aa





set.seed(123)

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)

library(Seurat)

library(dplyr)

library(Matrix)

library(methods)

library(RColorBrewer)

YDL#读入处理好的seurat对象
sweep.res.list <- paramSweep_v3(YDL, PCs = 1:20, sct = FALSE)
head(sweep.res.list)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
head(sweep.res.list)
bcmvn <- find.pK(sweep.stats)

annotations <- YDL@meta.data$seurat_clusters

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------

homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- YDL@meta.data$ClusteringResults

nExp_poi <- round(0.075*nrow(YDL@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))



## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

YDL <- doubletFinder_v3(YDL, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(YDL@meta.data)
YDL <- doubletFinder_v3(YDL, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_284", sct = FALSE)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")    #TSNE可视化
YDL_NEW <- subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$DF.classifications_0.25_0.09_284=="Singlet",]))
DimPlot(YDL_NEW,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")
DimPlot(YDL_NEW,reduction = "tsne",label = TRUE,pt.size = 1.5)

YDL_NEW_1 <- subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$DF.classifications_0.25_0.09_284=="Doublet",]))
DimPlot(YDL_NEW_1,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")
DimPlot(YDL_NEW_1,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by="DF.classifications_0.25_0.09_284")

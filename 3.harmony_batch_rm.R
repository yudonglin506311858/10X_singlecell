

set.seed(123)
library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(harmony)
#YDL<-readRDS("YDL.rds")
#肘部图（碎石图）,基于每个主成分对方差解释率的排名。
#看一个碎石图，看那个点主要在哪个地方没有什么下降趋势。然后选择这个数字作为PC数的值。
ElbowPlot(YDL)
YDL <- YDL %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(YDL, 'harmony')
harmony_embeddings[1:5, 1:5]
dim(harmony_embeddings)
head(harmony_embeddings)
p1 <- DimPlot(object = YDL, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = YDL, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
YDL <- YDL %>%
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  identity()


YDL <- YDL %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.4) %>%
  identity()






#绘制UMI的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)



library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","percent.mt"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = percent.mt))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



library("cowplot")
pdf("umi.pdf")
#绘制UMI的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
p2
plot_grid(p1,p2)
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)


#绘制基因的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
p2
plot_grid(p1,p2)
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)



#绘制线粒体的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","percent.mt"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = percent.mt))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
p2
plot_grid(p1,p2)
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","percent.mt"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = percent.mt))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)

dev.off()

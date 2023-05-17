
#Load demo dataset
library(CSOmapR)
library(CSOmapR.demo)
invisible(TPM)
invisible(LR)
invisible(labelData)


head(TPM[1:4,1:4])
head(LR)
head(labelData)

TPM=as.data.frame(YDL[["RNA"]]@counts) %>% apply(2,function(x){x/sum(x) * 10000})
#https://blog.csdn.net/m0_58549466/article/details/125730805


YDL[["celltype"]]<-Idents(YDL)
head(YDL@meta.data)
x = Idents(object = YDL)
head(x)
write.csv(x,"labelData.csv")
labelData<-read.csv("labelData.csv")
colnames(labelData)<-c("cells","labels")
# y<-YDL@meta.data$seurat_clusters
# head(y)
# x<-as.data.frame(x)
# head(x)
# y<-as.data.frame(y)
# head(y)
# new<-cbind(x,y)
# head(new)
# colnames(new)<-c("cells","labels")
# head(new)

#Calculate optimized 3D coordinates
affinityMat = getAffinityMat(TPM, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)
coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
#Visualization(by 3D density)
require(dplyr)
# arrange data
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)

head(density_obj)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
p_3Ddensity = plot3D(cellinfo_tbl, color_by = "labels", title = "3D density")
head(cellinfo_tbl$labels)
head(p_3Ddensity)
#Get significance
signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)
head(signif_results$qvalue)
head(contribution_list)
#When dealing with large dataset
#We provide two options: optimize coordinates through tSNE(BH algorithm), or downsample the original huge dataset first.


section_z0 = cellinfo_tbl[which(cellinfo_tbl$z < 1 & cellinfo_tbl$z > -1),]
ggplot(data = section_z0, aes(x = x, y = y, color = labels)) +
  geom_point()+theme_bw()
# under development
# Citation
# Ren, X., Zhong, G., Zhang, Q., Zhang, L., Sun, Y., and Zhang, Z. (2020). Reconstruction of cell spatial organization from single-cell RNA sequencing data based on ligand-receptor mediated self-assembly. Cell Res.
#https://github.com/lijxug/CSOmapR


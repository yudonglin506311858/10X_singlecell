#12:normal epithelial cells
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("EPCAM","KRT8","KRT18"),cols = c("gray", "red"))#actin
# fibroblasts (COL1A1+)
#4: liver tumor fibroblasts
#0:rectal cancer fibroblasts
#5:rectal cancer fibroblasts
#13:common tumor fibroblasts
#11:normal fibroblasts
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("COL1A1"),cols = c("gray", "red"))#actin
#:endothelial cells (CLDN5+)
#6,7,8,14
#2:rectal cancer endothelial cells
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("CLDN5"),cols = c("gray", "red"))#actin
#T cells (CD3D+)
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("CD3D"),cols = c("gray", "red"))#actin
#B cells (CD79A+)
#9
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("CD79A"),cols = c("gray", "red"))#actin
#myeloid cells (LYZ+)
#10
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("LYZ"),cols = c("gray", "red"))#actin


Idents(YDL) <- YDL@meta.data$seurat_clusters
current.cluster.ids <- c(0:12)
new.cluster.ids <- c(
  "fibroblasts",#0!
  "epithelial cells",#1!
  "endothelial cells",#2!
  "fibroblasts",#3!
  "endothelial cells",#4!
  "epithelial cells",#5!
  "myeloid cells",#6!
  "B cells",#7!
  "epithelial cells_normal",#8!
  "endothelial cells",#9!
  "unknown",#10
  "fibroblasts",#11!
  "fibroblasts"#12!
)
names(new.cluster.ids) <- levels(YDL)
YDL <- RenameIdents(YDL, new.cluster.ids)
pdf("细胞亚群注释.pdf",width = 6,height = 6)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "tsne", label = F, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = F, pt.size = 1.5)
dev.off()
pdf(file="umap__patient_rename.pdf",width=6.5,height=6)
UMAPPlot(object = YDL, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(object=YDL,label = TRUE,pt.size = 1.5, reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

pdf("细胞类型注释_patient1-B123.pdf")
#T细胞 CD2,CD3D,TRAC,TRBC2
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("CD2","CD3D","TRAC","TRBC2"),cols = c("gray", "red"))#actin
#Plasma细胞：JCHAIN,MZB1,CD79A,IGHG1
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("JCHAIN","MZB1","CD79A","IGHG1"),cols = c("gray", "red"))#actin
#Mononuclear phagocytes：LYZ,CD14,VCAN,C1QA,CD68,CDC
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("LYZ","CD14","VCAN","CD68","CDC"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("LYZ","CD14","VCAN","C1QA","CD68","CDC"),cols = c("gray", "red"))#actin
#Erythrocytes：HBB,HBA1,ALAS2,SNCA,CA1
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("HBB","HBA1","ALAS2","SNCA"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("HBB","HBA1","ALAS2","SNCA","CA1"),cols = c("gray", "red"))#actin
#Neutrophils：CSF3R,CXCR2,FCGR3B
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("CSF3R","CXCR2","FCGR3B"),cols = c("gray", "red"))#actin
#Fibroblasts：DCN,COL1A1,COL1A2
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("DCN","COL1A1","COL1A2"),cols = c("gray", "red"))#actin
#Epithelial cells：EPCAM,CLDN4,CLDN7,TFF3,REG4
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("EPCAM","CLDN4","CLDN7","TFF3","REG4"),cols = c("gray", "red"))#actin
#Granulocytes-monocyte progenitor cells：LYZ,MPO,AZU1,ELANE
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("LYZ","MPO","AZU1","ELANE"),cols = c("gray", "red"))#actin
#Platelets：PPBP,TUBB1,PF4
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("PPBP","TUBB1","PF4"),cols = c("gray", "red"))#actin
#Mast cells：TPSAB1,TPSB2,CPA3
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("TPSAB1","TPSB2","CPA3"),cols = c("gray", "red"))#actin
#Basophils：CLC,GATA2,CPA3
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("CLC","GATA2","CPA3"),cols = c("gray", "red"))#actin
dev.off()



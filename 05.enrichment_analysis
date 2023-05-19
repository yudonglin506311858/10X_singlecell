
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


library(dplyr)
library(clusterProfiler)
gid <- bitr(unique(YDL.markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(YDL.markers, gid, by=c('gene' = 'SYMBOL'))
x = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG')

write.csv(markers,file="marker-用于富集分析.csv")
write.csv(x,file="kegg_marker-细胞亚群.csv")
dotplot(x, showCategory =10)+ theme(axis.text.x = element_text(angle=45, hjust=1)) 
dotplot(x, showCategory =15)+ theme(axis.text.x = element_text(angle=45, hjust=1)) 


x = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichGO', OrgDb= 'org.Hs.eg.db')

write.csv(x,file="go_marker-细胞亚群.csv")
dotplot(x, showCategory =10)+ theme(axis.text.x = element_text(angle=45, hjust=1)) 
dotplot(x, showCategory =15)+ theme(axis.text.x = element_text(angle=45, hjust=1)) 
#https://mp.weixin.qq.com/s/nvOXM2hpZq3UI_Npjhk8HQ


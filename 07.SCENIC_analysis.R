#人类
#load in the motifannotation this will load it into your environment but as the name in which is given to the list argument
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

#rename the motif annnotion by attributing it to the variable that is in the error
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

set.seed(123)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
##==分析准备==##
dir.create("SCENIC")
dir.create("SCENIC/int")
#scRNA <-sce.all
scRNA <-readRDS("allsample_processed.rds")
#rm(YDL)
setwd("./SCENIC") 
##准备细胞meta信息
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
#colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA),1000)
scRNAsub <- scRNA[,subcell]
saveRDS(scRNAsub, "scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
# dim(exprMat)
# [1] 32285  1000
#head(cellInfo)
##设置分析环境
mydbDIR <- "/data/yudonglin/reference/scenic/"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=60,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC")
saveRDS(scenicOptions, "int/scenicOptions.rds")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间
##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
# scenicOptions <- initializeScenic(org="hgnc", 
#                                   nCores=9,
#                                   dbDir=mydbDIR, 
#                                   dbs = mydbs,
#                                   datasetTitle = "erythrogenesis")
##推断转录调控网络（regulon）
runSCENIC_2_createRegulons(scenicOptions)
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "erythrogenesis")
##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
saveRDS(exprMat_all,"exprMat_all.rds")
# rm(list=ls())
# scenicOptions<-readRDS("scenicOptions.rds")
# exprMat_all<-readRDS("exprMat_all.rds")
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

#runSCENIC_3_scoreCells(scenicOptions, exprMat=log2(as.matrix(scRNA@assays$RNA@counts)+1))



# #使用shiny互动调整阈值
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
# savedSelections <- shiny::runApp(aucellApp)
# #保存调整后的阈值
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)
savehistory("scenic_code.txt")





#小鼠

set.seed(123)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
##==分析准备==##
dir.create("SCENIC")
dir.create("SCENIC/int")
#scRNA <-sce.all
scRNA <-readRDS("allsample_processed.rds")
#rm(YDL)
setwd("./SCENIC") 
##准备细胞meta信息
cellInfo <- data.frame(scRNA@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
#colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
#为了节省计算资源，随机抽取1000个细胞的数据子集
subcell <- sample(colnames(scRNA),1000)
scRNAsub <- scRNA[,subcell]
saveRDS(scRNAsub, "scRNAsub.rds")
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
# dim(exprMat)
# [1] 32285  1000
#head(cellInfo)
##设置分析环境

#load in the motifannotation this will load it into your environment but as the name in which is given to the list argument
data(list="motifAnnotations_mgi_v9", package="RcisTarget")

#rename the motif annnotion by attributing it to the variable that is in the error
motifAnnotations_mgi <- motifAnnotations_mgi_v9

mydbDIR <- "/data/yudonglin/reference/scenic/"
mydbs <- c("mm9-500bp-upstream-7species.mc9nr.feather",
           "mm9-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="mgi", 
                                  nCores=60,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC")
saveRDS(scenicOptions, "int/scenicOptions.rds")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
#这一步消耗的计算资源非常大，个人电脑需要几个小时的运行时间
##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
#推断转录调控网络（regulon）
runSCENIC_2_createRegulons(scenicOptions)
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量
scenicOptions <- initializeScenic(org="mgi", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC")
##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
saveRDS(exprMat_all,"exprMat_all.rds")
# rm(list=ls())
# scenicOptions<-readRDS("scenicOptions.rds")
# exprMat_all<-readRDS("exprMat_all.rds")
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

#runSCENIC_3_scoreCells(scenicOptions, exprMat=log2(as.matrix(scRNA@assays$RNA@counts)+1))



# #使用shiny互动调整阈值
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
# savedSelections <- shiny::runApp(aucellApp)
# #保存调整后的阈值
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)
savehistory("scenic_code.txt")



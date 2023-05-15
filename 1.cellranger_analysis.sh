

下载软件
wget -O cellranger-7.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1660521741&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NjA1MjE3NDF9fX1dfQ__&Signature=I1ClmYpGC95JhgfEilOLlMOtkbvNUm3uVXqjnmHEz2vFYPVhMF1aZIk-p8NHV37YUABtlrfX5wYUygwMtYMQ6ZpC3bBNiY29rEWf0QQCUsHCxh23v6-q2TQMYlckMntqIeCD63CkxjWLdX85nEgepLpmLKDSH2gSIMxP0PyhAA1TzHdVDIkPcGa2b7ShosdFK7CQXvXJ~nNDYa4AIrNq9Ku1C1RerSNjEPYHd2Bkz7b1pSoEM69rutHgIA43fkCtuW9xX3hXt-lSYtqBNRYyLABqJelx3JQ4lWVufUdEh-YYJZuA2u0RtRF4uTe96KdENkafebKuByH3O9PaJV9Hag__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
 
下载官方基因组：
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz
 

软件安装：
tar -xzvf cellranger-7.0.0.tar.gz
export PATH=/data/yudonglin/software/cellranger-7.0.0:$PATH
 
#cellranger testrun --id=tiny



cd /data3/yudonglin/ren

Individual Accession	Individual Name	Gender	Sample Accession	Sample Name	Sample Description
HRI015279	INDIVIDUAL_HRA000063_15279	male	HRS015279	20170608M	10X of human: liver metastasized tumor from rectal cancer
HRI015280	INDIVIDUAL_HRA000063_15280	male	HRS015280	20170608P	10X of human: primary rectal cancer
HRI015281	INDIVIDUAL_HRA000063_15281	male	HRS015281	20171109N	10X of human: adjacent non-tumor tissue
HRI015282	INDIVIDUAL_HRA000063_15282	male	HRS015282	20171109T	10X of human: liver tumor



mv HRR016237_f1.fastq.gz HRR016237_S1_L001_R1_001.fastq.gz
mv HRR016237_r2.fastq.gz HRR016237_S1_L001_R2_001.fastq.gz

#配置环境，成功安装：
export PATH=/data/yudonglin/software/cellranger-7.0.0:$PATH
#开始运行：
cellranger count --id=HRR016237 \
                   --transcriptome=/data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A\
                   --fastqs=/data3/yudonglin/ren/ \
                   --sample=HRR016237 \
                   --localcores=40 \
                   --localmem=64



export OPENBLAS_NUM_THREADS=1 
conda activate velocity 
velocyto run10x -m /data/yudonglin/reference/singcell/velocity/hg19_rmsk.gtf /data3/yudonglin/ren/HRR016237 /data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A/genes/genes.gtf 
 


cd /data3/yudonglin/ren
mv HRR016238_f1.fastq.gz HRR016238_S1_L001_R1_001.fastq.gz
mv HRR016238_r2.fastq.gz HRR016238_S1_L001_R2_001.fastq.gz

#配置环境，成功安装：
export PATH=/data/yudonglin/software/cellranger-7.0.0:$PATH
#开始运行：
cellranger count --id=HRR016238 \
                   --transcriptome=/data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A\
                   --fastqs=/data3/yudonglin/ren/ \
                   --sample=HRR016238 \
                   --localcores=40 \
                   --localmem=64



export OPENBLAS_NUM_THREADS=1 
conda activate velocity 
velocyto run10x -m /data/yudonglin/reference/singcell/velocity/hg19_rmsk.gtf /data3/yudonglin/ren/HRR016238 /data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A/genes/genes.gtf 
 



cd /data3/yudonglin/ren
mv HRR016239_f1.fastq.gz HRR016239_S1_L001_R1_001.fastq.gz
mv HRR016239_r2.fastq.gz HRR016239_S1_L001_R2_001.fastq.gz

#配置环境，成功安装：
export PATH=/data/yudonglin/software/cellranger-7.0.0:$PATH
#开始运行：
cellranger count --id=HRR016239 \
                   --transcriptome=/data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A\
                   --fastqs=/data3/yudonglin/ren/ \
                   --sample=HRR016239 \
                   --localcores=40 \
                   --localmem=64


export OPENBLAS_NUM_THREADS=1 
conda activate velocity 
velocyto run10x -m /data/yudonglin/reference/singcell/velocity/hg19_rmsk.gtf /data3/yudonglin/ren/HRR016239 /data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A/genes/genes.gtf 
 



cd /data3/yudonglin/ren
mv HRR016240_f1.fastq.gz HRR016240_S1_L001_R1_001.fastq.gz
mv HRR016240_r2.fastq.gz HRR016240_S1_L001_R2_001.fastq.gz

#配置环境，成功安装：
export PATH=/data/yudonglin/software/cellranger-7.0.0:$PATH
#开始运行：
cellranger count --id=HRR016240 \
                   --transcriptome=/data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A\
                   --fastqs=/data3/yudonglin/ren/ \
                   --sample=HRR016240 \
                   --localcores=40 \
                   --localmem=64



export OPENBLAS_NUM_THREADS=1 
conda activate velocity 
velocyto run10x -m /data/yudonglin/reference/singcell/velocity/hg19_rmsk.gtf /data3/yudonglin/ren/HRR016240 /data/yudonglin/reference/singcell/refdata-gex-GRCh38-2020-A/genes/genes.gtf 
 


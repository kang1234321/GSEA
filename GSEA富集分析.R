############################### GSEA 富集分析 ##################################
#将得到的差异基因（包含logFC）进行富集分析
# 加载GSEA分析所需要的包
#https://zhuanlan.zhihu.com/p/667669372
setwd('E:\\生信\\结肠癌转移\\单细胞\\GSE178318')
library(clusterProfiler) # GSEA富集/基因集读取
library(dplyr)
library(tidyverse)
# 导入差异分析后的数据，以便后续使用logFC进行基因排序
#load("./gsea_data/DEG_limma.Rdata")
rm(list = ls())
DEG_limma <- read.csv('上皮细胞亚型差异1.csv')
row.names(DEG_limma) <- DEG_limma$X
colnames(DEG_limma)[3] <- 'logFC'
head(DEG_limma)

DEG_limma<- DEG_limma %>%
  filter(str_detect(cluster, "Stem-like cell"))#提取数据框cluster列包含Stem-like cell的行
# 加载基因集，基因集介绍往下滑！
geneSet_h <- read.gmt("h.all.v2023.2.Hs.symbols.gmt")#使用GSEA网站下载的hallmark基因集
head(geneSet_h)
# 接下来我们进行基因排序
geneList <- DEG_limma$logFC                 # 获取GeneList
names(geneList) <- DEG_limma$X     # 对GeneList命名
head(geneList)
class(geneList)
geneList <- sort(geneList, decreasing = T)  # 从高到低排序
# 开始GSEA富集分析
set.seed(123)
GSEA_enrichment <- GSEA(geneList,                 # 排序后的gene
                        TERM2GENE = geneSet_h,    # 基因集
                        pvalueCutoff = 1,      # P值阈值
                        minGSSize = 10,           # 最小基因数量
                        maxGSSize = 1000,         # 最大基因数量
                        eps = 0,                  # P值边界
                        pAdjustMethod = "BH"
                      )     # 校正P值的计算方法


result <- data.frame(GSEA_enrichment)
#subset_data <- result [result$pvalue < 0.05, ]#筛选p<0.05的行
dim(GSEA_enrichment@result)
write.csv(result,'stemlike亚型H富集分析.csv',row.names = T)
# 展示富集到的通路，我们这里选择展示前15个
dotplot(GSEA_enrichment, showCategory = 10, color = "pvalue")



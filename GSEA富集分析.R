##############家人们投一下币，点一下关注吧！！
##############QQ群：1121444037  为了防止有人进来打广告，请不要使用小号加群！！！管理员可能不会给你通过！
##############由于QQ群开放讨论，请注意防骗。
##############2023-1-5
##############所有课程都是紧密联系的，如果没看过之前的课程请不要观看这节课！！！！
##############以后的课用到的R语言技能都是基础篇和中级篇学过的！没有任何R基础纯粹想套代码的话，建议放弃。
##############出课程的目的一是给大家快速入门R语言，二用来纪录我自己的学习内容。
##############------------------------------分割线----------------------------------------------############
##############--------------------------项目实战课7-GSEA富集分析-----------------------------########################################
##############GSEA代码地址：https://www.jianshu.com/p/bb7442bb8cad
##############              https://zhuanlan.zhihu.com/p/358168557#
##############设置工作目录
setwd('E:\\R study\\Ｒ语言学习分享 (1)\\P20_GSEA富集分析')
##############
library(clusterProfiler)
library(dplyr)
##############
deg.data <- data.table::fread('DEG.csv',data.table = F)#获取差异基因表达谱
###############获取需要富集分析的基因
#############修改列名
colnames(deg.data)[3] <- 'logFC'
colnames(deg.data)[7] <- 'adj.P.Val'
############新增一列，区分上调和下调的基因
deg.data$Group <- 'Not-significant'
deg.data$Group[which((deg.data$adj.P.Val) < 0.0001 & (deg.data$logFC > 5) )] <- 'Up'#设置筛选条件
deg.data$Group[which((deg.data$adj.P.Val) < 0.0001 & (deg.data$logFC < -5) )] <- 'Down'

############获取基因列表和LogFC值
deg <- deg.data %>% 
                    filter(Group %in% c('Down','Up')) %>% 
                    dplyr::select(1,3)#选中deg.data文件的1  3列

############Symbol ID转换为ENTREZID id ---原先表达矩阵差异分析的基因时sympleID，则对下载的基因集用转换，若是 ENTREZID id 基因集文件则不需要转换          
gene <-  bitr(deg$id, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")  
############基因去重复
gene <- distinct(gene,SYMBOL,.keep_all = T)#.对gene文件的symble列操作，keep_all = T为保留后面的列

############获取LogFC添加大gene中，使用merge函数
colnames(gene)[1] <- 'id'  #数据框列名改成相同
gene <- deg %>% #提供数据的文件
                merge(gene,by='id')#交叉匹配
  
############获取genelist
genelist <- gene$logFC
############使用entrezid ID命名或者symbol ID
names(genelist) <- gene$id#此处使用的是sympleID
############对genelist排序
genelist <- sort(genelist,decreasing = T)#从大到小排列
genelist 

test <- c(7,1,9,3,8)
a <- sort(test) ## 默认升序
a <- sort(test,decreasing = T)#decreasing = T为降序排列
a
############读入需要富集分析的GMT文件。
############网页地址：http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
############以KEGG为例
kegmt <- read.gmt('c2.cp.kegg.v2022.1.Hs.symbols.gmt')#读取gmt格式文件，此为下载的基因集 ，第一列为基因通路名，第二列为网址或空，第三列及之后为基因集中基因

############正式进行GSEA分析
kkgesa <- GSEA(genelist,TERM2GENE =kegmt)#使用GSEA函数，对genelist文件中的基因进行以kegmt基因集的富集分析

############
#ES：富集分数
#NES：标准化后的富集分数
#NOM p-val：是对 ES 的统计学分析，用来表征富集结果的可信度
#FDR q-val：是多重假设检验校正之后的 p-value，即对NES可能存在的假阳性结果的概率估计。

############查看结果
result <- as.data.frame(kkgesa)#结果转化为数据框

############绘图
dotplot(kkgesa)#绘制气泡图
###########以P值出图
dotplot(kkgesa,color='pvalue')
###########显示通路激活或者抑制
library(ggplot2)
dotplot(kkgesa,split='.sign')+facet_grid(~.sign)

###########特定通路作图
###########enrichplot https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
###########clusterplofier  

library(enrichplot)
result$ID#可查看基因富集的通路
gseaplot2(kkgesa,"KEGG_RETINOL_METABOLISM",color = "red",pvalue_table = T)
gseaplot2(kkgesa,"KEGG_DRUG_METABOLISM_CYTOCHROME_P450",color = "red",pvalue_table = T)
gseaplot2(kkgesa,7,color = "red",pvalue_table = T)#数字7代表将第7条通路进行富集
gseaplot2(kkgesa, geneSetID = 1:4)#表示将第1-4条通路富集到一张图中


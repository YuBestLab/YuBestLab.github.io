#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

allFile="all.txt"      #所有基因的差异结果文件
gmtFile="c6.all.v2023.1.Hs.symbols.gmt"      #基因集文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\19model_differ\\04功能性分析\\04原癌基因信号通路\\TCGA")     #设置工作目录

#读取输入文件,并对输入文件进行整理
rt=read.table(allFile, header=T, sep="\t", check.names=F)
rt=rt[order(rt[,"logFC"],decreasing=T),]
logFC=as.vector(rt[,"logFC"])
names(logFC)=as.vector(rt[,1])

#读入基因集文件
gmt=read.gmt(gmtFile)

#对排序好的基因进行GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

#输出ICD高表达组富集的图形
kkUp=kkTab[kkTab$NES>0,]
termNum=5     #设置展示通路的数目，展示前5个富集最显著的通路
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]      #获取展示通路的名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in ICD high group")
  pdf(file="GSEA.highrisk.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
}

#输出ICD低表达组富集的图形
kkDown=kkTab[kkTab$NES<0,]
termNum=5     #设置展示通路的数目，展示前5个富集最显著的通路
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]      #获取展示通路的名称
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in ICD low group")
  pdf(file="GSEA.lowrisk.pdf", width=7, height=5.5)
  print(gseaplot)
  dev.off()
}


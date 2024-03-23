#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #引用包
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\05_1TMB\\02GSH代谢基因的瀑布图")      #设置工作目录

#读取突变基因文件
geneRT=read.table("BP_gene.txt", header=T, sep="\t", check.names=F)
gene=geneRT[,1]

#绘制瀑布图
pdf(file="oncoplot.pdf", width=8, height=7.5)
maf=read.maf(maf="input.maf")
oncoplot(maf=maf, genes=gene, fontSize=0.5, draw_titv=T)
dev.off()

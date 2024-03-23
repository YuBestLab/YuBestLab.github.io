

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)

expFile="GSHexp.txt"         #表达数据文件
clusterFile="cluster.txt"     #分型结果文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\02_2PCA")      #设置工作目录
#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))

#读取分组的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample,drop=F]
data=t(data)
cluster=cluster[sameSample,,drop=F]
#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab.xls", quote=F, sep="\t")

#读取分型文件
PRGcluster=as.vector(cluster[,1])

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
prgCluCol=bioCol[1:length(levels(factor(PRGcluster)))]

#可视化
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], PRGcluster=PRGcluster)
PCA.mean=aggregate(PCA[,1:2], list(PRGcluster=PCA$PRGcluster), mean)
pdf(file="PCA.pdf", width=6.5, height=5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = PRGcluster)) +
  scale_colour_manual(name="PRGcluster", values =prgCluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$PRGcluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



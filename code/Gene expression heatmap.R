#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)

expFile="GSHexp.txt"          #表达数据文件
ClusterFile="cluster.txt"     #分型的结果文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\03GSHgroup\\TCGA")     #设置工作目录

#读取表达数据文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=log2(data+1)

#读取分型的结果文件
Cluster=read.table(ClusterFile, header=T, sep="\t", check.names=F, row.names=1)
Type=Cluster[order(Cluster$Cluster),,drop=F]
exp=t(data[row.names(Type),])

#设置热图注释的颜色
ann_colors=list()
bioCol=c("lightskyblue2", "lightcoral")
names(bioCol)=c("C1", "C2")
ann_colors[["Cluster"]]=bioCol

#绘制热图
pdf(file="heatmap.pdf", width=7, height=5.5)
pheatmap(exp,
         annotation=Type,
         annotation_colors=ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=8,
         fontsize_row=8,
         fontsize_col=8)
dev.off()

#对样品进行分组
median1=median(exp[,Type$Cluster=="C1"])     #获取分型1表达的中位值
median2=median(exp[,Type$Cluster=="C2"])     #获取分型2表达的中位值
if(median1>median2){
  Type$Group=ifelse(Type$Cluster=="C1", "GSH high", "GSH low")
}else{
  Type$Group=ifelse(Type$Cluster=="C1", "GSH low", "GSH high")
}
#输出分组的结果
Type=rbind(ID=colnames(Type), Type)
write.table(Type, file="GSHgroup.txt", sep="\t", quote=F, col.names=F)

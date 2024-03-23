#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("pheatmap")


#引用包
library(limma)
library(ggplot2)
library(pheatmap)

conGroup="GSH low"         #对照组
treatGroup="GSH high"      #实验组
expFile="LGG.txt"       #表达数据文件
cluFile="GSHgroup.txt"     #分组的结果文件
logFCfilter=1              #logFC过滤条件(logFC=0.585,差异倍数1.5倍;logFC=1,差异2倍;logFC=2,差异4倍)
fdrFilter=0.05             #fdr过滤条件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\05groupdiff")     #设置工作目录

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
Type=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(Type))
data=data[,sameSample,drop=F]
Type=Type[sameSample,,drop=F]

#提取不同分组的样品
low=Type[Type$Group==conGroup,,drop=F]          #提取ICD低表达组的样品
high=Type[Type$Group==treatGroup,,drop=F]       #提取ICD高表达组的样品
dataLow=data[,row.names(low)]
dataHigh=data[,row.names(high)]
data=cbind(dataLow, dataHigh)
data=data[rowMeans(data)>0.1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pvalue=pvalue))
  }
}
pvalue=outTab[,"pvalue"]
fdr=p.adjust(as.numeric(as.vector(pvalue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)
write.table(outTab, file="all.txt", sep="\t", row.names=F, quote=F)

#输出差异分析的结果
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="diff.txt", sep="\t", row.names=F, quote=F)

#输出差异基因的表达数据
diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]), data[as.vector(outDiff[,1]),])
write.table(diffExp, file="diffGeneExp.txt", sep="\t", col.names=F, quote=F)

#绘制差异基因的热图
geneNum=50    #定义展示基因的数目
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep(conGroup,conNum),rep(treatGroup,treatNum))
Type=factor(Type, levels=c(conGroup, treatGroup))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

#绘制差异基因的火山图
library(dplyr)
library(ggplot2)   
#install.packages("gt")
library(gt)
library(ggrepel)


setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\05groupdiff\\01提取分组间的差异基因") 
logFCfilter=1              #logFC过滤条件(logFC=0.585,差异倍数1.5倍;logFC=1,差异2倍;logFC=2,差异4倍)
fdrFilter=0.05             #fdr过滤条件
rt<-read.table("all.txt", header = T, check.names = F, row.names = 1)
pdf(file="vol1.pdf", width=6, height=5)
#定义显著性
rt$Significant=ifelse((rt$fdr<fdrFilter&abs(rt$logFC)>1), ifelse(rt$logFC>1,"Up","Down"), "Not")
###标记前10的基因
#top_10 <-bind_rows(
 # rt %>%
 #   filter(Significant == 'Up') %>%
  #  arrange(fdr, desc(abs(logFC))) %>%
   # head(5),
 # rt %>%
  #  filter(Significant == 'Down') %>%
   # arrange(fdr, desc(abs(logFC))) %>%
    #head(5))
#top_20 %>% gt()

#绘制火山图
P1 <- ggplot(rt, aes(logFC, -log10(fdr)))+
  geom_hline(yintercept = -log10(0.05), colour="grey", linetype="dashed") +
  geom_vline(xintercept = c(-1,1), colour="grey", linetype="dashed") +
  geom_point(aes(col=Significant), size = 1,alpha = 0.3)+
  scale_color_manual(values=c("mediumturquoise", "grey", "indianred1"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  theme_bw()
#P2 <- P1 +geom_text_repel(data = top_10,
                          #aes(logFC, -log10(fdr), label = rownames(top_10)),box.padding=unit(0.2, "lines"),
                          #force = 1,segment.color = "dimgrey", point.padding=unit(0.2, "lines"),
                          #size=1, segment.size=0.1,segment.alpha=0.5)
#出图
print(P1)
dev.off()

####
#绘制火山图
pdf(file="vol3.pdf", width=6, height=5)
P1 <- ggplot(rt, aes(logFC, -log10(fdr)))+
  geom_hline(yintercept = -log10(0.05), colour="black", linetype="dashed", size = 0.4) +
  geom_vline(xintercept = c(-1,1), colour="black", linetype="dashed", size = 0.4) +
  geom_point(aes(col=Significant), size = 0.8)+
  scale_color_manual(values=c("cornflowerblue", "grey", "tomato"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  theme_bw()

#出图
print(P1)
dev.off()




#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("pheatmap")


#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

expFile="LGG.txt"               #表达数据文件
clusterFile="GSHgroup.txt"      #分型结果文件
gmtFile="c2.cp.kegg.v2023.1.Hs.symbols.gmt"                   #基因集文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\05groupdiff\\03分组间GSVA分析")      #设置工作目录

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

#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#读取cluster文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
#Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
#gsvaCluster=cbind(gsvaCluster, Project)

#差异分析
adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$Cluster)
comp=combn(levels(factor(allType)), 2)
for(i in 1:ncol(comp)){
  #样品分组
  treat=gsvaCluster[gsvaCluster$Cluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$Cluster==comp[1,i],]
  data=rbind(con, treat)
  #对通路进行差异分析
  Type=as.vector(data$Cluster)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有通路的差异情况
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出显著的差异
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  
  #设置热图注释的颜色
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(allType))
  ann_colors[["cluster"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
  
  #绘制差异通路热图
  termNum=20     #设置显示通路的数目
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".heatmap.pdf"), width=10, height=6)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
           cluster_cols =F,
           show_colnames = F,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize = 8,
           fontsize_row=6,
           fontsize_col=8)
  dev.off()
}



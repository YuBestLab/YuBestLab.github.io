#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")
#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("ComplexHeatmap")


#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library("ggpubr")
library(ComplexHeatmap)

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=0.05      #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\05groupdiff\\02差异基因的GO, KEGG分析")       #设置工作目录
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#提取差异基因的名词, 将基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
##
####各别基因对应两个entrzID,因此导出基因手段改为后重新读入
write.csv(gene, "gene.csv")####各别基因对应两个entrzID,因此导出基因手段改为后重新读入
gen<-read.csv("gene.csv", row.names = 1, header = T)
gene=as.character(gen[,1])
###
#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#定义显示GO的数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="barplot.pdf", width=12, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=60, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=12, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=60, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

###########绘制GO圈图###########
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
data=GO[order(GO$p.adjust),]
datasig=data[data$p.adjust<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图的主体部分
pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制显著性pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


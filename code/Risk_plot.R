#install.packages("pheatmap")


library(pheatmap)         #引用包
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\12riskplot\\01无多因素")      #设置工作目录

#定义风险曲线的函数
bioRiskPlot=function(inputFile=null, project=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
  rt=rt[order(rt$riskScore),]      #根据病人风险得分对样品进行排序
 
  #绘制风险曲线
  riskClass=rt[,"Risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=paste0(project, ".riskScore1.pdf"), width=7, height=4)
  plot(line, type="p", pch=19, cex=0.7,
       xlab="Patients (increasing risk socre)",
       ylab="Risk score",
       col=c(rep("mediumturquoise",lowLength),rep("indianred1",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk","Low risk"),bty="n",pch=19,col=c("indianred1","mediumturquoise"),cex=0.7)
  dev.off()
  
  #绘制生存状态图
  color=as.vector(rt$fustat)
  color[color==1]="indianred1"
  color[color==0]="mediumturquoise"
  pdf(file=paste0(project, ".survStat1.pdf"), width=7, height=4)
  plot(rt$futime, pch=19, cex=0.7,
       xlab="Patients (increasing risk socre)",
       ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("indianred1","mediumturquoise"),cex=0.7)
  abline(v=lowLength,lty=2)
  dev.off()
  
  #定义热图注释的颜色
  ann_colors=list()
  bioCol=c("mediumturquoise", "indianred1")
  names(bioCol)=c("low", "high")
  ann_colors[["Risk"]]=bioCol
  
  #绘制风险热图
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(Risk=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=paste0(project, ".heatmap1.pdf"), width=7, height=4)
  pheatmap(rt1, 
           annotation=annotation,
           annotation_colors = ann_colors, 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = F,
           color = colorRampPalette(c(rep("mediumturquoise",3.5), "white", rep("indianred1",3.5)))(50),
           scale="row",
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}

#调用函数，绘制风险曲线
bioRiskPlot(inputFile="risk_CCGA693.txt", project="CCGA")
bioRiskPlot(inputFile="risk.TCGA.txt", project="TCGA")


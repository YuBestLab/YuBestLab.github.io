inputFile="Gene_CNV_matrix.txt"     #输入文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\05_2基因拷贝数变异\\02棒棒糖图")      #设置工作目录

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
GAIN=rowSums(rt> 0)         #获取拷贝数增加的样品数目
LOSS=rowSums(rt< 0)         #获取拷贝数缺失的样品数目
GAIN=GAIN/ncol(rt)*100      #计算拷贝数增加的百分率
LOSS=LOSS/ncol(rt)*100      #计算拷贝数缺失的百分率
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#绘制图形
data.max = apply(data, 1, max)
pdf(file="CNVfreq1.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1, cex=0.7)
dev.off()

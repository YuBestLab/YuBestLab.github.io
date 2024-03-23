#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)

GroupFile="GSHgroup.txt"     #分组的结果文件
cliFile="time.txt"           #生存数据文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\04gene Cluster survival\\TCGA\\OS")      #设置工作目录

#读取输入文件
Group=read.table(GroupFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(Group), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], Group[sameSample,,drop=F])
rt[,"Group"]=factor(rt[,"Group"], levels=c("GSH low","GSH high"))

#比较ICD高低表达组的生存差异
length=length(levels(factor(rt$Group)))
diff=survdiff(Surv(futime, fustat) ~ Group, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Group, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("mediumturquoise","indianred1","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,   # 增加置信区间
                   pval=pValue,
                   pval.size=6,
                   legend.title="Group",
                   legend.labs=levels(factor(rt[,"Group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by=2,
                   palette=bioCol,
                   surv.median.line="hv",    # 添加中位生存时间线
                   risk.table=TRUE,          # 添加风险表
                   fontsize=3,                         # 指定风险表和累积事件表的字体大小
                   risk.table.title="",
                   ggtheme = theme_bw(),
                   risk.table.height=.25)

#输出图形
pdf(file="survival.pdf", width=6.5, height=5, onefile=FALSE)
print(surPlot)
dev.off()

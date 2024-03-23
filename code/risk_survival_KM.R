#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\11model_survival\\无多因素_旧")     #设置工作目录

#定义生存分析的函数
bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组的生存差异，得到差异显著性的pvalue
  diff=survdiff(Surv(futime, fustat) ~ Risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     ylab="Overall survival",
                     break.time.by = 2,
                     palette=c("indianred1","mediumturquoise"),
                     risk.table=TRUE,          # 添加风险表
                     fontsize=3,                         # 指定风险表和累积事件表的字体大小
                     risk.table.title="",
                     ggtheme = theme_bw(),
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, width=6, height=5, onefile=FALSE)
  print(surPlot)
  dev.off()
}

#调用函数,绘制生存曲线
bioSurvival(inputFile="risk.TCGA.txt", outFile="survival.TCGA1.pdf")
####bioSurvival(inputFile="risk_geo.txt", outFile="risk_CCGA693.txt")

####CCGA693的图
rt=read.table("risk_CCGA693.txt", header=T, sep="\t", check.names=F, row.names = 1)
#比较高低风险组的生存差异，得到差异显著性的pvalue
diff=survdiff(Surv(futime, fustat) ~ Risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt)

#绘制生存曲线
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Risk",
                   legend.labs=c("High risk", "Low risk"),
                   xlab="Time(years)",
                   ylab="Overall survival",
                   break.time.by = 2,
                   palette=c("indianred1","mediumturquoise"),
                   risk.table=TRUE,          # 添加风险表
                   fontsize=3,                         # 指定风险表和累积事件表的字体大小
                   risk.table.title="",
                   ggtheme = theme_bw(),
                   risk.table.height=.25)
#输出图形
pdf(file="CCGA1_693.pdf", width=6, height=5, onefile=FALSE)
print(surPlot)
dev.off()

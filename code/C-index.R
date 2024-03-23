#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("survcomp")

#install.packages("survival")
#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(survival)
library(survcomp)
library(ggplot2)
library(ggpubr)

inputFile="risk.models.txt"      #输入文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\10model\\04与其它模型的比较\\02Cindex")       #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)


##########绘制C-index图形##########
#获取每个模型的C-index值
df=data.frame()
for(i in colnames(rt)[3:ncol(rt)]){
  cindex=concordance.index(x=rt[,i], surv.time=rt$futime, surv.event=rt$fustat,method="noether")
  df=rbind(df, cbind(i,sprintf("%.03f",cindex$c.index)))
}
colnames(df)=c("signature", "cindex")
df[,"cindex"]=as.numeric(df[,"cindex"])

#绘制C-index的柱状图
color=rainbow(nrow(df),alpha=0.75)
p=ggbarplot(df, x="signature", y="cindex", fill="signature",
            xlab="", ylab="C-index", add = "none",
            palette=color,
            label=T, legend="")
p=p+rotate_x_text(50)
p=p+ylim(0,round(max(df[,"cindex"])+0.15,1))
#输出图形
pdf(file="C-index.pdf", width=6, height=5)
print(p)
dev.off()


##########绘制RMS图形##########
outdata=list()
legendsname=c()
#获取每个模型的RMS值
for(i in 3:ncol(rt)){
  #数据整理
  OS=Surv(rt$futime, rt$fustat)
  marker=rt[,i]
  marker.pp<-seq(from=0, to=1, length=100)
  marker.qq<-quantile(marker,marker.pp)
  fitdat.df<-data.frame(marker=marker)
  newdat.df<-data.frame(marker=marker.qq)
  
  #计算RMS
  cox.model<-coxph(OS~marker, data=fitdat.df)
  rms.calc <-summary(survfit(cox.model, newdata=newdat.df))
  rms.mean <-rms.calc$table[,"rmean"]
  name=colnames(rt)[i]
  HR=sprintf("%.03f", summary(cox.model)$conf.int[,"exp(coef)"])
  HR.95L=sprintf("%.03f", summary(cox.model)$conf.int[,"lower .95"])
  HR.95H=sprintf("%.03f", summary(cox.model)$conf.int[,"upper .95"])
  pvalue=summary(cox.model)$coefficients[,"Pr(>|z|)"]
  p=ifelse(pvalue<0.001,"p<0.001",paste0("p=",sprintf("%.03f",pvalue)))
  
  #获取图例名称
  legendsname=c(legendsname,paste0(name,", HR:",HR,"(",HR.95L,"-",HR.95H,"), ",p))     #图例
  outdata[[name]]= data.frame(marker.pp,rms.mean)
}

#绘制RMS图形
alldata=do.call("rbind",outdata)
xlim2=max(alldata$rms.mean)
pdf(file="RMS.pdf", width=6, height=6)
par(las=1)
plot(1,xlim=c(0,1),ylim=c(0,xlim2),type="n",xlab="Percentile of scores",ylab="RMS (years)")
names=names(outdata)
for(i in 1:length(outdata)){
  namei=names[i]
  outdatai=outdata[[namei]]
  points(outdatai$marker.pp,outdatai$rms.mean,col=color[i],pch=20,cex=0.8)
}
legend("bottomleft",legend=legendsname,col=color,pch=20,bty="n",cex=1)
dev.off()

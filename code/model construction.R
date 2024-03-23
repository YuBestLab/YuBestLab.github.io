#install.packages("glmnet")
#install.packages("survival")


#引用包
library("glmnet")
library("survival")

trainFile="TCGA.uniSigExp.txt"    #单因素显著基因的表达文件
testFile="CCGA693_expTime.txt"        #GEO数据库输入文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\10model")     #设置工作目录

#读取train组数据文件
rt=read.table(trainFile, header=T, sep="\t", check.names=F, row.names=1)
rt$futime[rt$futime<=0]=0.003

#构建lasso回归模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
#绘制lasso回归的图形
pdf("lasso.lambda.pdf")
plot(fit, xvar="lambda", label=TRUE)
dev.off()
#绘制交叉验证的图形
cvfit=cv.glmnet(x, y, family="cox", maxit=1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#找到交叉验证误差最小的点，并且输出模型公式
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

#输出TCGA数据库的风险文件
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
Risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.TCGA.txt",sep="\t",quote=F,row.names=F)

#输出GEO数据库的风险文件
rt=read.table(testFile, header=T, sep="\t", row.names=1)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
Risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="risk.GEO.txt",sep="\t",quote=F,row.names=F)


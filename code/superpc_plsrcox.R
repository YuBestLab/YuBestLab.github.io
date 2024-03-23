# 加载需要使用的R包
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(plsRcox)
library(superpc)
source(file.path("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\21十种机器学习组合\\单个学习方法作图\\Codes", "ML.R"))

##################################
#### 准备工作 ####
##################################

# 加载数据集
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\21十种机器学习组合\\单个学习方法作图")
tcga <- read.table("TCGA.uniSigExp.txt", header = T,sep = "\t", quote = "", check.names = FALSE, row.names = 1)
#GSE57303 <- read.table("C:\\Users\\lenovo\\Desktop\\全文复现\\GSE57303.txt", header = T, sep = "\t", quote = "", check.names = F)
#GSE62254 <- read.table("C:\\Users\\lenovo\\Desktop\\全文复现\\GSE62254.txt", header = T, sep = "\t", quote = "", check.names = F)

# 生成包含三个数据集的列表
mm <- list(TCGA = tcga)

# 数据标准化
mm <- lapply(mm,function(x){
  x[,-c(1:2)] <- scaleData(data = tcga[,-c(1:2)], centerFlags = T, scaleFlags = T)
  return(x)})

# TCGA作为训练集
est_data <- mm$TCGA
# GEO作为验证集
#val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:2)]
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]
#val_dd_list <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', pre_var)]})


rf_nodesize <- 5
seed <- 123


##################################
#### 5.plsRcox####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$OS.time,status=est_dd$OS),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
par(mfrow=c(1,2))
DR_coxph(est_dd$OS.time, est_dd$OS,plot = T)
DR_coxph(est_dd$OS.time, est_dd$OS, scaleY = FALSE , plot = T)
##################################
#### 6.superpc####
##################################

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=3, 
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

superpc.plotcv(cv.fit)

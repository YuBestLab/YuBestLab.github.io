######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager::install(c("limma", "car", "ridge", "preprocessCore", "genefilter", "sva", "biomaRt"))
#BiocManager::install("GenomicFeatures", "maftools", " stringr", "org.Hs.eg.db")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

install.packages("oncoPredict")

#引用包
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="LGG.txt"     #表达数据文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\17药物敏感性分析\\oncoPredict")     #设置工作目录

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

#读取药物敏感性文件
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#药物敏感性
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = data,
              batchCorrect = 'eb',      #"eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056




#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")


#引用包
library(TwoSampleMR)
library(ieugwasr)
inputFile="exposureID_gsh.txt"      #暴露数据id文件
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\22孟德尔随机化\\01提取暴露因素")     #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F,)

#对暴露数据ID进行循环
outTab=data.frame()
for(id in rt$ID) {
	expoData=extract_instruments(id,
                                 p1 = 5e-08, p2 = 5e-07,
                                 clump = T,
                                 kb = 10000, r2 = 0.01)
    outTab=rbind(outTab, expoData)
}
write.csv(outTab, file="exposure_data_8_8_10k_0.01.csv", row.names=F)


#install.packages("ggpubr")

library(scales)
library(ggpubr)              
inputFile="geneCoef.txt"        
outFile="geneCoef.pdf"      
setwd("D:\\生信分析课题\\泛癌数据\\谷胱甘肽途径\\10model\\01model\\coef作图")    

#读取
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
x=colnames(rt)[1]
y=colnames(rt)[2]
#colnames(rt)=c("Name","Value")
rt$Regulate=factor(ifelse(rt$Coef<0, "Down", "Up"), levels = c("Up", "Down"))

#绘制
pdf(file=outFile,width=8,height=5)
ggbarplot(rt, x="Gene", y="Coef", fill = "Regulate",  
          color = "black",
          palette = c("indianred1","mediumturquoise"), 
          sort.val = "desc", sort.by.groups = FALSE, 
          x.text.angle=90, 
          xlab = x, ylab = y, 
          legend.title="Regulate", rotate=TRUE, ggtheme = theme_minimal())+ #rotatex设置x/y对调
          theme_bw()+
          scale_y_continuous(breaks = breaks_width(0.2, offset = -0.2), labels = label_number(accuracy = 0.1))+ #自定义坐标轴刻度，以0.2为间隔
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


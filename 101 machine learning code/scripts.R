
# 本代码适用于：
# 存在较多预后意义不明的变量，需同时进行变量筛选，模型筛选和预后模型构建的场景
# 
# 代码输出结果包括：
# 101种算法在训练集上获得的模型(model.rds)，在测试集(或包括训练集)上的评估结果(Cindex_mat)，以及评估结果的热图展示(Cindex.pdf)。

# 设置工作路径
work.path <- "D:\\生信分析课题\\泛癌数据\\SKCM免疫疗法\\03_多种机器学习组合\\0代码示例\\01机器学习\\方法三_取GEO交集预后基因"; setwd(work.path) 

# 设置其他路径
code.path <- file.path(work.path, "Codes") # 存放脚本
data.path <- file.path(work.path, "InputData") # 存在输入数据（需用户修改）
res.path <- file.path(work.path, "Results") # 存放输出结果
fig.path <- file.path(work.path, "Figures") # 存放输出图片

# 如不存在这些路径则创建路径
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# install.packages("randomForestSRC")
# install.packages("snowfall")

# 加载需要使用的R包
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))

# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiCox表示抽取其他模型所用到的变量并建立多变量cox模型
FinalModel <- c("panML", "multiCox")[2]

## Training Cohort ---------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，表达谱需有一定变异性，以免建模过程报错）
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 训练集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 测试集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因

# 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
## data: 需要表达谱数据（行为样本，列为基因） 
## cohort：样本所属队列，为向量，不输入值时默认全表达矩阵来自同一队列
## centerFlag/scaleFlags：是否将基因均值/标准差标准化为1；
##        默认参数为NULL，表示不进行标准化；
##        为T/F时，表示对所有队列都进行/不进行标准化
##        输入由T/F组成的向量时，按顺序对队列进行处理，向量长度应与队列数一样
##        如centerFlags = c(F, F, F, T, T)，表示对第4、5个队列进行标准化，此时flag顺序应当与队列顺序一致
##        如centerFlags = c("A" = F, "C" = T, "B" = F)，表示对队列C进行标准化，此时不要求flag顺序与data一致
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# lapply(split(as.data.frame(Test_set), Test_surv$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 5 # 筛选变量数目的最小阈值
timeVar = "OS.time"; statusVar = "OS" # 定义需要考虑的结局事件，必须出现在Train_surv以及Test_surv中

## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

set.seed(seed = 123) # 设置建模种子，使得结果可重复
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 机器学习方法
                                 Train_expr = Train_set, # 训练集有潜在预测价值的变量
                                 Train_surv = Train_surv, # 训练集生存数据
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的生存变量，必须出现在Train_surv中
}
preTrain.var[["simple"]] <- colnames(Train_expr)

model <- list() # 初始化模型结果列表
set.seed(seed = 123) # 设置建模种子，使得结果可重复
for (method in methods){ # 循环每一种方法组合
  # method <- "CoxBoost+plsRcox" # [举例]若遇到报错，请勿直接重头运行，可给method赋值为当前报错的算法来debug
  cat(match(method, methods), ":", method, "\n") # 输出当前方法
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method)
  
  selected.var = preTrain.var[[method[1]]]
  # 如果筛选出的变量小于阈值，则该算法组合无意义，置空（尤其针对以RSF筛选变量的情况，需在ML脚本中尝试调参）
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method[2], # 用于构建最终模型的机器学习方法
                                  Train_expr = Train_expr[, selected.var], # 训练集有潜在预测价值的变量
                                  Train_surv = Train_surv, # 训练集生存数据
                                  mode = "Model",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                  classVar = classVar)  # 用于训练的生存变量，必须出现在Train_surv中
  }
  
  # 如果最终筛选出的变量小于阈值，则该算法组合也无意义，置空
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
saveRDS(model, file.path(res.path, "model.rds")) # 保存所有模型输出

# 当要求最终模型为多变量cox时，对模型进行更新
if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量cox模型，从而以cox回归系数和特征表达计算单样本风险得分
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
    return(tmp)
  })
}
saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) # 保存最终以多变量cox拟合所筛选变量的模型

## Evaluate the model -----------------------------------------------------

# 读取已保存的模型列表（请根据需要调整）
model <- readRDS(file.path(res.path, "model.rds")) # 若希望使用各自模型的线性组合函数计算得分，请运行此行
# model <- readRDS(file.path(res.path, "coxmodel.rds")) # 若希望使用多变量cox模型计算得分，请运行此行

methodsValid <- names(model) # 取出有效的模型（变量数目小于阈值的模型视为无效）

# 根据给定表达量计算样本风险评分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = rbind.data.frame(Train_set,Test_set), # 4.0更新
                                    type = "lp") # 同原文，使用linear Predictor计算得分
  
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  # 数据框有两列，包含算法以及算法所筛选出的变量
write.table(fea_df, file.path(res.path, "fea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = Test_set, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_surv = Test_surv, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                  Train_expr = Train_set, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_surv = Train_surv, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                  Train_name = "TCGA", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  #Train_expr = NULL,
                                  #Train_surv = NULL, 
                                  cohortVar = "Cohort", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  timeVar = timeVar, # 用于评估的生存时间，必须出现在Test_surv中；这里是OS.time
                                  statusVar = statusVar) # 用于评估的生存状态，必须出现在Test_surv中；这里是OS
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth = 1.5; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 3, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())

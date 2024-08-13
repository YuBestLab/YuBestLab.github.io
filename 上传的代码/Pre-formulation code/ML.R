# all ML algorithms
RunML <- function(method, Train_expr, Train_surv, mode = "Model", timeVar = "OS.time", statusVar = "OS", ...){
  # for example: Enet [alpha=0.4]
  method = gsub(" ", "", method) # 去除参数中的空格，得到Enet [alpha=0.4]
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # get name of ML algorithm, e.g., Enet
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # get parameter of ML algorithm, e.g., alpha=0.4
  
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "StepCox" = list("direction" = method_param),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_expr), " Variables")
  
  args = list("Train_expr" = Train_expr,
              "Train_surv" = Train_surv,
              "mode" = mode,
              "timeVar" = timeVar, "statusVar" = statusVar)
  args = c(args, method_param)
  
  obj <- do.call(what = paste0("Run", method_name),
                 args = args) 
  
  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_expr, Train_surv, mode, timeVar, statusVar, alpha){
  cv.fit = cv.glmnet(x = Train_expr,
                     y = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]),
                     family = "cox", alpha = alpha, nfolds = 10)
  fit = glmnet(x = Train_expr,
               y = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]),
               family = "cox", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  RunEnet(Train_expr, Train_surv, mode, timeVar, statusVar, alpha = 1)
}

RunRidge <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  RunEnet(Train_expr, Train_surv, mode, timeVar, statusVar, alpha = 0)
}

RunStepCox <- function(Train_expr, Train_surv, mode, timeVar, statusVar, direction){
  fit <- step(coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                    data = as.data.frame(Train_expr)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunsurvivalSVM <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  fit = survivalsvm(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                    data= as.data.frame(Train_expr),
                    gamma.mu = 1, opt.meth = "ipop")
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunCoxBoost <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  pen <- optimCoxBoostPenalty(time = Train_surv[[timeVar]],
                              status = Train_surv[[statusVar]],
                              x = Train_expr,
                              trace = F, start.penalty=500, parallel = F)
  cv.res <- cv.CoxBoost(time = Train_surv[[timeVar]],
                        status = Train_surv[[statusVar]],
                        x = Train_expr,
                        maxstepno=500, K=10, type="verweij", penalty=pen$penalty)
  fit <- CoxBoost(time = Train_surv[[timeVar]],
                  status = Train_surv[[statusVar]],
                  x = Train_expr,
                  stepno = cv.res$optimal.step,
                  penalty = pen$penalty)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSuperPC <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  data <- list(x = t(scale(Train_expr)),
               y = Train_surv[[timeVar]],
               censoring.status = Train_surv[[statusVar]],
               featurenames = colnames(Train_expr))
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) # default
  cv.fit <- suppressWarnings(superpc.cv(fit, data,
                                        n.threshold = 20, # default
                                        n.fold = 10,
                                        n.components = 3,
                                        min.features = 5,
                                        max.features = nrow(data$x),
                                        compute.fullcv = TRUE,
                                        compute.preval = TRUE))
  fit$threshold <- cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])]
  fit$data <- data
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRcox <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  data <- list(x = Train_expr,
               time = Train_surv[[timeVar]],
               status = Train_surv[[statusVar]])
  cv.plsRcox.res = cv.plsRcox(data = data,
                              nt=10, verbose = FALSE)
  fit <- plsRcox(Xplan = data$x,
                 time = data$time,
                 event = data$status,
                 nt = as.numeric(cv.plsRcox.res[5]),
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRSF <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  rf_nodesize = 5 # 如果RSF无有效的变量筛选结果，则需要反复调参
  fit <- rfsrc(formula = formula(paste0("Surv(", timeVar, ", ", statusVar, ")", "~.")),
               data = cbind(Train_expr, Train_surv[,c(timeVar, statusVar)]), # 3.1更新
               ntree = 1000, nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunGBM <- function(Train_expr, Train_surv, mode, timeVar, statusVar){
  fit <- gbm(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
             data = as.data.frame(Train_expr),
             distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
             data = as.data.frame(Train_expr),
             distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = colnames(Train_expr)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "coxnet" = rownames(coef(fit))[which(coef(fit)[, 1] != 0)], # 本身没有筛选变量的功能，但是可以舍去模型中系数为0的变量
    "coxph" = names(coef(fit)), # 逐步回归可以对变量进行筛选
    "survivalsvm" = fit$var.names, # survivalSVM对变量没有筛选功能，所以默认使用所有变量
    "CoxBoost" = names(coef(fit)[abs(coef(fit)) > 0]), # CoxBoost同样不具备筛选变量的功能，因此舍去模型中系数为0的变量
    "superpc" = names(fit$feature.scores)[abs(fit$feature.scores) > 0.5],# superpct同样不具备筛选变量的功能，因此舍去模型中特征得分小于0.5的特征
    "plsRcoxmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0], # plsRcoxmodel同样不具备筛选变量的功能，因此舍去模型中系数为0的变量
    "rfsrc" = var.select(fit, verbose = F)$topvars, # rfsrc可以对变量进行筛选
    "gbm" = rownames(summary.gbm(fit, plotit = FALSE))[summary.gbm(fit, plotit = FALSE)$rel.inf > 0] # gbm通过舍去重要性为0的变量来进行变量筛选
  ))
  return(Feature)
}

CalRiskScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, fit$subFeature]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "coxnet"      = predict(fit, type = 'link', as.matrix(new_data)),
    "coxph"       = predict(fit, type = 'lp', as.data.frame(new_data)),
    "survivalsvm" = predict(fit, as.data.frame(new_data))$predicted,
    "CoxBoost"    = predict(fit, type = "lp", as.data.frame(new_data)),
    "superpc"     = superpc.predict(object = fit,
                                    data = fit$data,
                                    newdata = list(x = t(scale(as.matrix(new_data)))),
                                    threshold = fit$threshold,
                                    n.components = 1)$v.pred,
    "plsRcoxmodel" = predict(fit, type = "lp", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted,
    "gbm"          = predict(fit, type = 'link', as.data.frame(new_data))
  ))
  RS = as.vector(RS)
  names(RS) = rownames(new_data)
  return(RS)
}

RunEval <- function(fit, 
                    Test_expr = NULL, 
                    Test_surv = NULL, 
                    Train_expr = NULL, 
                    Train_surv = NULL, 
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    timeVar = "OS.time", 
                    statusVar = "OS"){
  
  if(!is.element(cohortVar, colnames(Test_surv))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  } 
  
  if((!is.null(Train_expr)) & (!is.null(Train_surv))) {
    new_data <- rbind.data.frame(Train_expr[, fit$subFeature],
                                 Test_expr[, fit$subFeature])
    
    if(!is.null(Train_name)) {
      Train_surv$Cohort <- Train_name
    } else {
      Train_surv$Cohort <- "Training"
    }
    colnames(Train_surv)[ncol(Train_surv)] <- cohortVar
    Test_surv <- rbind.data.frame(Train_surv[,c(cohortVar, timeVar, statusVar)],
                                  Test_surv[,c(cohortVar, timeVar, statusVar)])
    Test_surv[,1] <- factor(Test_surv[,1], 
                            levels = c(unique(Train_surv[,cohortVar]), setdiff(unique(Test_surv[,cohortVar]),unique(Train_surv[,cohortVar]))))
  } else {
    new_data <- Test_expr[, fit$subFeature]
  }
  
  RS <- CalRiskScore(fit = fit, new_data = new_data)
  
  Predict.out <- Test_surv
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  f <- as.formula(paste0("Surv(", timeVar,",",statusVar,")~RS"))
  unlist(lapply(Predict.out, function(data){
    unname(summary(coxph(formula = f,
                         data = data))$concordance["C"])
  }))
}

SimpleHeatmap <- function(Cindex_mat = NULL, 
                          avg_Cindex = NULL, 
                          CohortCol = NULL, 
                          barCol = NULL,
                          col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  Heatmap(as.matrix(Cindex_mat), name = "C-index",
          right_annotation = row_ha, 
          top_annotation = col_ha,
          col = col, 
          rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 不进行聚类，无意义
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}

# 3.1更新
standarize.fun <- function(indata, centerFlag, scaleFlag) {  
  scale(indata, center=centerFlag, scale=scaleFlag)
}

# 4.0更新
scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }
  
  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  # lapply(out.data, function(x) summary(apply(x, 2, var)))
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}
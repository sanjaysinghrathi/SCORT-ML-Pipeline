#### 1 #### Validation function
expression_data <- exprs(public_myeset)
metadata <- public_pData[,1:3]
val_cohort_size <- "Balanced Cohort"
combat_type <- "response+tn"
validate_GSE87211<-function(expression_data, metadata, val_cohort_size , combat_type){
  library(Biobase)
  ## Load data
  load(paste0(path_to_discovery_data,"DiscoveryBig.RData"))
  expression_data <- as.matrix(expression_data)
  val1_testing_fix_myeset <- ExpressionSet(assayData = expression_data)
  val1_testing_fix_pData <- metadata
  colnames(val1_testing_fix_pData) <- c("response","Pretreatment.T.Stage","Pretreatment.N.Stage")
  phenoData(val1_testing_fix_myeset) <- new("AnnotatedDataFrame", data = val1_testing_fix_pData)
  row.names(val1_testing_fix_pData)==colnames(val1_testing_fix_myeset)
  val1_testing_fix_e <- exprs(val1_testing_fix_myeset)
  
  ######### Get common gene data for training and testing sets for validation 1
  val1_training_myeset <- raw_training_myeset
  val1_training_pData <- raw_training_pData
  val1_training_e <- exprs((val1_training_myeset))
  dim(val1_training_e)
  dim(val1_testing_fix_e)
  common <- intersect(row.names(deg_list), row.names(val1_testing_fix_e))
  #common <- intersect(row.names(val1_training_e), row.names(val1_testing_fix_e))
  length(common)
  val1_training_e <- val1_training_e[common,]
  val1_testing_fix_e <- as.matrix(val1_testing_fix_e[common,])
  colnames(val1_testing_fix_e) <- row.names(val1_testing_fix_pData)
  dim(val1_training_e)
  dim(val1_testing_fix_e)
  val1_training_myeset <- ExpressionSet(assayData = val1_training_e)
  phenoData(val1_training_myeset) <- new("AnnotatedDataFrame", data = val1_training_pData)
  val1_testing_fix_myeset <- ExpressionSet(assayData = val1_testing_fix_e)
  phenoData(val1_testing_fix_myeset) <- new("AnnotatedDataFrame", data = val1_testing_fix_pData)
  
  ## train data
  aab <- data.frame(val1_training_pData$response)
  colnames(aab) <- c("response")
  pp <- val1_training_e[intersect(row.names(deg_list),row.names(val1_training_e)),]
  training <- cbind(t(as.data.frame(pp)), aab)
  rm(aab)
  rm(pp)
  training$response = factor(training$response)
  dim(training)
  training$response = factor(training$response, levels = c("responder","nonresponder"))
  
  ########################################### feature selection using random forest
  #load(paste0(path_to_validation_data,"public-data-gene.RData"))
  edata <- exprs(public_myeset)
  if(length(setdiff(row.names(edata),row.names(expression_data)))==0){
    ## prepare data
    n <- dim(training)[2]
    aab <- data.frame(training[,n])
    colnames(aab) <- c("response")
    training1 <- cbind(training[,bestest1], aab)
    rm(aab)
    rm(n)
  }else{
    set.seed(123)
    boruta.train <- Boruta(response~., data = training, doTrace = 2,  maxRuns= 3000, pValue = 0.05)
    print(boruta.train)
    final.boruta <- TentativeRoughFix(boruta.train)
    print(final.boruta)
    rm(bestest1)
    set.seed(123)
    bestest1 <- data.frame(significant_probes =getSelectedAttributes(final.boruta, withTentative = TRUE))
    rm(boruta.train)
    rm(final.boruta)
    
    library(tidyr)
    bestest1 <- separate_rows(bestest1, 1, sep = "`")
    bestest1 <- bestest1[bestest1$significant_probes!="",]
    bestest1 <- bestest1[[1]]
    #final_deg_list <- final_deg_list[as.character(bestest1),]
    
    ## prepare data
    n <- dim(training)[2]
    aab <- data.frame(training[,n])
    colnames(aab) <- c("response")
    training1 <- cbind(training[,bestest1], aab)
    rm(aab)
    rm(n)
  }
  
  ###########################################################    Feature based Gradient boosting algo
  # Tune random forest
  set.seed(567)
  numFolds <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,  allowParallel = TRUE, verbose=TRUE , search = "grid")
  grid <- expand.grid(interaction.depth = c(1, 3, 5, 7, 9), 
                      n.trees = (1:30)*50, 
                      shrinkage = c(0.1,0.01, 1, 10, 100),
                      n.minobsinnode = c(10,20))
  
  set.seed(567)
  FGBM <- train(response ~ ., training1, method='gbm', metric="Accuracy", trControl = numFolds, tuneGrid=grid)
  #models[[12]] <- FGBM
  
  ################ Prepare test data for validation (remove batch effect)
  val1_testing_e <- val1_testing_fix_e
  val1_testing_myeset <- val1_testing_fix_myeset
  val1_testing_pData <- val1_testing_fix_pData
  
  sample_names <- row.names(val1_testing_pData)
  val1_testing_pData$responsenstage <- paste(val1_testing_pData$Pretreatment.T.Stage ,val1_testing_pData$Pretreatment.N.Stage, sep=".")
  val1_testing_pData <- val1_testing_pData %>% 
    mutate(responsenstage = recode(responsenstage, 
                                   "1.0" = "Block1",
                                   "1.1" = "Block2",
                                   "1.2" = "Block2",
                                   "2.0" = "Block1",
                                   "2.1" = "Block2",
                                   "2.2" = "Block2",
                                   "3.0" = "Block3",
                                   "3.1" = "Block4",
                                   "3.2" = "Block4",
                                   "4.0" = "Block3",
                                   "4.1" = "Block4",
                                   "4.2" = "Block4"))
  val1_testing_pData$responsetstage <- val1_testing_pData$response
  val1_testing_pData <- val1_testing_pData %>% 
    mutate(responsetstage = recode(responsetstage, 
                                   "nonresponder" = "0",
                                   "responder" = "1"))
  
  val1_testing_pData$responsetstage <- paste(val1_testing_pData$responsetstage ,val1_testing_pData$responsenstage, sep="")
  val1_testing_pData <- as.data.frame(lapply(val1_testing_pData, function (x) if (is.factor(x)) factor(x) else x))
  #val1_testing_pData$responsetstage <- as.numeric(as.character(val1_testing_pData$responsetstage))
  phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
  
  val1_testing_pData$binarystage <- paste(val1_testing_pData$Pretreatment.T.Stage ,val1_testing_pData$Pretreatment.N.Stage, sep=".")
  val1_testing_pData <- val1_testing_pData %>% 
    mutate(binarystage = recode(binarystage, 
                                "1.0" = "T12-N0",
                                "1.1" = "T12-N12",
                                "1.2" = "T12-N12",
                                "2.0" = "T12-N0",
                                "2.1" = "T12-N12",
                                "2.2" = "T12-N12",
                                "3.0" = "T34-N0",
                                "3.1" = "T34-N12",
                                "3.2" = "T34-N12",
                                "4.0" = "T34-N0",
                                "4.1" = "T34-N12",
                                "4.2" = "T34-N12"))
  
  val1_testing_pData$responsebinarystage <- paste(val1_testing_pData$response ,val1_testing_pData$binarystage, sep=".")
  val1_testing_pData <- as.data.frame(lapply(val1_testing_pData, function (x) if (is.factor(x)) factor(x) else x))
  #val1_testing_pData$responsetstage <- as.numeric(as.character(val1_testing_pData$responsetstage))
  phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
  
  val1_testing_pData$continousstage <- paste(val1_testing_pData$Pretreatment.T.Stage ,val1_testing_pData$Pretreatment.N.Stage, sep=".")
  val1_testing_pData <- val1_testing_pData %>% 
    mutate(continousstage = recode(continousstage, 
                                   "1.0" = "1",
                                   "1.1" = "2",
                                   "1.2" = "2",
                                   "2.0" = "1",
                                   "2.1" = "2",
                                   "2.2" = "2",
                                   "3.0" = "3",
                                   "3.1" = "4",
                                   "3.2" = "4",
                                   "4.0" = "3",
                                   "4.1" = "4",
                                   "4.2" = "4"))
  
  val1_testing_pData$responsecontinousstage <- paste(val1_testing_pData$response ,val1_testing_pData$continousstage, sep=".")
  val1_testing_pData <- as.data.frame(lapply(val1_testing_pData, function (x) if (is.factor(x)) factor(x) else x))
  #val1_testing_pData <- val1_testing_pData[, c(1:3,6:11,4:5)]
  row.names(val1_testing_pData) <- sample_names
  #val1_testing_pData <- cbind(val1_testing_pData, grampian_pData[,5:8])
  phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
  val1_testing_e <- exprs(val1_testing_myeset)
  
  if(length(setdiff(row.names(edata),row.names(expression_data)))==0){
    combat_type <- "response+tn"
  }else{
    combat_type <- "tn"
  }
  combat <- combat_type
  if(combat=="response+tn"){
    if (!requireNamespace("sva", quietly = TRUE)){
      BiocManager::install("sva")}
    
    library(sva)
    full_em <- cbind(val1_training_e, val1_testing_e)
    full_pData <- rbind(val1_training_pData[,1:5], val1_testing_pData[,c(1:5)])
    batch <- c(rep(1, dim(val1_training_e)[2]),rep(2, dim(val1_testing_e)[2]))
    mod <- model.matrix(~as.factor(response), data=full_pData)
    set.seed(567)
    combat_edata3 <- ComBat(dat=full_em, batch=batch, mod = mod, par.prior = TRUE, ref.batch=1)
    #training_em_new <- combat_edata3[,1:dim(training_em)[2]]
    val1_testing_e <- as.matrix(combat_edata3[, (dim(val1_training_e)[2]+1):dim(full_em)[2]])
    colnames(val1_testing_e) <- row.names(val1_testing_pData)
    val1_testing_myeset <- ExpressionSet(assayData = val1_testing_e)
    phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
    ## training data her
    val1_training_e <- combat_edata3[, 1:(dim(val1_training_e)[2])]
    val1_training_myeset <- ExpressionSet(assayData = val1_training_e)
    phenoData(val1_training_myeset) <- new("AnnotatedDataFrame", data = val1_training_pData)
  }else{
    if (!requireNamespace("sva", quietly = TRUE)){
      BiocManager::install("sva")}
    
    library(sva)
    full_em <- cbind(val1_training_e, val1_testing_e)
    full_pData <- rbind(val1_training_pData[,1:5], val1_testing_pData[,c(1:5)])
    batch <- c(rep(1, dim(val1_training_e)[2]),rep(2, dim(val1_testing_e)[2]))
    mod <- model.matrix(~as.factor(response), data=full_pData)
    set.seed(567)
    combat_edata3 <- ComBat(dat=full_em, batch=batch, par.prior = TRUE, ref.batch=1)
    #training_em_new <- combat_edata3[,1:dim(training_em)[2]]
    val1_testing_e <- as.matrix(combat_edata3[, (dim(val1_training_e)[2]+1):dim(full_em)[2]])
    colnames(val1_testing_e) <- row.names(val1_testing_pData)
    val1_testing_myeset <- ExpressionSet(assayData = val1_testing_e)
    phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
    ## training data her
    val1_training_e <- combat_edata3[, 1:(dim(val1_training_e)[2])]
    val1_training_myeset <- ExpressionSet(assayData = val1_training_e)
    phenoData(val1_training_myeset) <- new("AnnotatedDataFrame", data = val1_training_pData)
  }
  
  
  #### create 10 iterations of validation balanced set
  val1_testing_set <- list() 
  if (!requireNamespace("raster", quietly = TRUE)){
    BiocManager::install("raster")}
  set.seed(45)
  seeds <- raster::sampleInt(1000, 10, replace = FALSE)
  for(iter in 1:10){
    rand <- seeds[iter]
    val1_testing_myeset_temp <- val1_testing_myeset
    val1_testing_pData_temp <- val1_testing_pData
    val1_testing_e_temp <- val1_testing_e
    testsampling <- val_cohort_size
    if(testsampling =="Balanced Cohort"){
      ## downsampling of training data
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block1")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block1")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat1 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat1 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block2")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block2")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat2 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat2 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block3")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block3")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat3 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat3 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block4")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block4")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat4 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat4 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      
      bal_e7 <- val1_testing_e_temp
      add_e7 <- val1_testing_e_temp[, c(cat1, cat2, cat3,cat4)]
      dim(add_e7)
      add_e7[1:5,1:5]
      dim(add_e7)
      add_e7[1:5,1:5]
      
      bal_pData <- val1_testing_pData_temp
      add_pData <- val1_testing_pData_temp[c(cat1, cat2, cat3,cat4),]
      dim(add_pData)
      add_pData[1:5, 1:5]
      
      val1_testing_pData_temp <- add_pData
      val1_testing_myeset_temp <- ExpressionSet(assayData = add_e7)
      phenoData(val1_testing_myeset_temp) <- new("AnnotatedDataFrame", data = val1_testing_pData_temp)
      val1_testing_e_temp <- exprs(val1_testing_myeset_temp)
      rm(bal_e7)
      rm(bal_pData)
      rm(add_e7)
      rm(add_pData)
      rm(cat1)
      rm(cat2)
      rm(cat3)
      rm(cat4)
      rm(cat1a)
      rm(cat1b)
    }
    ## Prepare test data for validation 1
    ## test data
    aab <- data.frame(val1_testing_pData_temp$response)
    colnames(aab) <- c("response")
    testing <- cbind(t(as.data.frame(val1_testing_e_temp)), aab)
    rm(aab)
    testing$response = factor(testing$response)
    dim(testing)
    testing$response = factor(testing$response, levels = c("responder","nonresponder"))
    val1_testing_set[[iter]] <- testing
  }
  length(intersect(row.names(val1_testing_set[[1]]), row.names(val1_testing_set[[10]])))
  
  ## check accuracies for all validation set testing combinations
  exo <-0
  acc_scores <- c()
  for(iter in 1:10){
    testing <- val1_testing_set[[iter]]
    #pro_temp <- pro+iter
    exo_temp <- exo+iter
    model <- FGBM
    bestest1 <- data.frame(significant_probes =model$coefnames)
    bestest1 <- separate_rows(bestest1, 1, sep = "`")
    bestest1 <- bestest1[bestest1$significant_probes!="",]
    testing_new <- testing[,c(bestest1$significant_probes, "response")]
    training$response = factor(training$response, levels = c("responder", "nonresponder"))
    testing_new$response = factor(testing_new$response, levels = c("responder", "nonresponder"))
    
    predictions <- model %>% predict(testing_new)
    predictions
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_new$response))){
      ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
    }else if(nlevels(as.factor(predictions))>=nlevels(as.factor(testing_new$response))){
      levels(testing_new$response) <- levels(predictions)
      ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
    }else if(nlevels(as.factor(predictions))<=nlevels(as.factor(testing_new$response))){
      levels(predictions) <- levels(testing_new$response)
      ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
    }
    acc_scores[exo_temp] <- ext_accuracy$overall[1]
  }
  
  ## Get final Model
  model <- FGBM
  testing <- val1_testing_set[[which(acc_scores== sort.default(acc_scores)[6])[1]]]
  
  bestest1 <- data.frame(significant_probes =model$coefnames)
  bestest1 <- separate_rows(bestest1, 1, sep = "`")
  bestest1 <- bestest1[bestest1$significant_probes!="",]
  testing_new <- testing[,c(bestest1$significant_probes, "response")]
  val1_final_deg_list <- deg_list[bestest1$significant_probes,]
  
  int_accuracy <- max(model$results$Accuracy, na.rm = TRUE)
  predictions <- predict(model, newdata = testing_new, type = "prob")
  row.names(predictions) <- row.names(testing_new)
  predictions_prob <- data.frame(predictions)
  predictions <- predict(model, newdata = testing_new)
  predictions
  predictions_prob$Predicted_class <- predictions
  predictions_prob$Actual_class <- testing_new$response
  
  predictions <- model %>% predict(testing_new)
  predictions
  
  if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_new$response))){
    ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
  }else if(nlevels(as.factor(predictions))>=nlevels(as.factor(testing_new$response))){
    levels(testing_new$response) <- levels(predictions)
    ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
  }else if(nlevels(as.factor(predictions))<=nlevels(as.factor(testing_new$response))){
    levels(predictions) <- levels(testing_new$response)
    ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
  }
  
  ## Print External Accuracy
  ext_accuracy
  
  ## plotROC
  if(!requireNamespace("plotROC", quietly = TRUE)){
    BiocManager::install("plotROC")}
  library(plotROC)
  predictions_prob <- predictions_prob
  #predictions_prob <- predictions_prob[order(predictions_prob$responder),]
  predictions_prob$label <- predictions_prob$Actual_class
  sname <- row.names(predictions_prob)
  predictions_prob <- predictions_prob %>% 
    mutate(label = recode(label, 
                          "nonresponder" = 0,
                          "responder" = 1))
  row.names(predictions_prob) <- sname
  rocplot <- ggplot(predictions_prob, aes(m = responder, d = label))+ geom_roc(n.cuts=20,labels=FALSE)
  #rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") 
  roc_plot <- rocplot + 
    style_roc(theme = theme_grey) +
    theme(axis.text = element_text(colour = "blue")) +
    ggtitle("ROC Curve for Grampian") + 
    annotate("text", x = .55, y = .35, 
             label = paste("AUC =", round(calc_auc(rocplot)$AUC, 2))) +
    annotate("text", x = .55, y = .30, 
             label = paste("Cutoff Point =", 0.50)) +
    scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1))
  
  predictions_prob1 <- predictions_prob
  thresh <- 0.50
  predictions_prob1$Predicted_class1 <- factor(ifelse(predictions_prob1$responder > thresh, "responder", "nonresponder"), levels = c("responder", "nonresponder"))
  confusionMatrix(table(predictions_prob1$Predicted_class1, predictions_prob$Actual_class))
  
  SampleNumber <- c(1:dim(testing)[1])
  prob_plot <- predictions_prob %>% 
    ggplot(aes(x=SampleNumber, y=responder, col=Actual_class)) + 
    scale_color_manual(values=c("red", "black")) + 
    geom_point() + 
    geom_rug() + 
    ggtitle("Response scores for Testing set")+
    geom_hline(yintercept=thresh, linetype="dashed", color = "red")
  
  #rm(predictions_prob1)
  rm(thresh)
  rm(predictions_prob)
  
  ## clear
  rm(val1_testing_myeset_temp)
  rm(val1_testing_e_temp)
  rm(val1_testing_fix_pData)
  rm(full_em)
  rm(full_pData)
  rm(mod)
  rm(batch)
  rm(combat_edata3)
  rm(common)
  rm(combat)
  rm(testsampling)
  #
  return(list(prob_plot, roc_plot,ext_accuracy, predictions_prob1))
}




















































##### 2 Predictiom funtiuon
expression_data <- exprs(public_myeset)
metadata <- public_pData[,2:3]
val_cohort_size <- "Balanced Cohort"
combat_type <- "response+tn"
predict_GSE87211<- function(expression_data, metadata, val_cohort_size , combat_type){
  ## Load data
  load(paste0(path_to_discovery_data,"DiscoveryBig.RData"))
  expression_data <- as.matrix(expression_data)
  val1_testing_fix_myeset <- ExpressionSet(assayData = expression_data)
  val1_testing_fix_pData <- metadata
  colnames(val1_testing_fix_pData) <- c("Pretreatment.T.Stage","Pretreatment.N.Stage")
  phenoData(val1_testing_fix_myeset) <- new("AnnotatedDataFrame", data = val1_testing_fix_pData)
  row.names(val1_testing_fix_pData)==colnames(val1_testing_fix_myeset)
  val1_testing_fix_e <- exprs(val1_testing_fix_myeset)
  
  ######### Get common gene data for training and testing sets for validation 1
  val1_training_myeset <- raw_training_myeset
  val1_training_pData <- raw_training_pData
  val1_training_e <- exprs((val1_training_myeset))
  dim(val1_training_e)
  dim(val1_testing_fix_e)
  common <- intersect(row.names(deg_list), row.names(val1_testing_fix_e))
  #common <- intersect(row.names(val1_training_e), row.names(val1_testing_fix_e))
  length(common)
  val1_training_e <- val1_training_e[common,]
  val1_testing_fix_e <- as.matrix(val1_testing_fix_e[common,])
  colnames(val1_testing_fix_e) <- row.names(val1_testing_fix_pData)
  dim(val1_training_e)
  dim(val1_testing_fix_e)
  val1_training_myeset <- ExpressionSet(assayData = val1_training_e)
  phenoData(val1_training_myeset) <- new("AnnotatedDataFrame", data = val1_training_pData)
  val1_testing_fix_myeset <- ExpressionSet(assayData = val1_testing_fix_e)
  phenoData(val1_testing_fix_myeset) <- new("AnnotatedDataFrame", data = val1_testing_fix_pData)
  
  ## train data
  aab <- data.frame(val1_training_pData$response)
  colnames(aab) <- c("response")
  pp <- val1_training_e[intersect(row.names(deg_list),row.names(val1_training_e)),]
  training <- cbind(t(as.data.frame(pp)), aab)
  rm(aab)
  rm(pp)
  training$response = factor(training$response)
  dim(training)
  training$response = factor(training$response, levels = c("responder","nonresponder"))
  
  ########################################### feature selection using random forest
  #load(paste0(path_to_validation_data,"public-data-gene.RData"))
  edata <- exprs(public_myeset)
  if(length(setdiff(row.names(edata),row.names(expression_data)))==0){
    ## prepare data
    n <- dim(training)[2]
    aab <- data.frame(training[,n])
    colnames(aab) <- c("response")
    training1 <- cbind(training[,bestest1], aab)
    rm(aab)
    rm(n)
  }else{
    set.seed(123)
    boruta.train <- Boruta(response~., data = training, doTrace = 2,  maxRuns= 3000, pValue = 0.05)
    print(boruta.train)
    final.boruta <- TentativeRoughFix(boruta.train)
    print(final.boruta)
    rm(bestest1)
    set.seed(123)
    bestest1 <- data.frame(significant_probes =getSelectedAttributes(final.boruta, withTentative = TRUE))
    rm(boruta.train)
    rm(final.boruta)
    
    library(tidyr)
    bestest1 <- separate_rows(bestest1, 1, sep = "`")
    bestest1 <- bestest1[bestest1$significant_probes!="",]
    bestest1 <- bestest1[[1]]
    #final_deg_list <- final_deg_list[as.character(bestest1),]
    
    ## prepare data
    n <- dim(training)[2]
    aab <- data.frame(training[,n])
    colnames(aab) <- c("response")
    training1 <- cbind(training[,bestest1], aab)
    rm(aab)
    rm(n)
  }
  
  ###########################################################    Feature based Gradient boosting algo
  # Tune random forest
  set.seed(567)
  numFolds <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,  allowParallel = TRUE, verbose=TRUE , search = "grid")
  grid <- expand.grid(interaction.depth = c(1, 3, 5, 7, 9), 
                      n.trees = (1:30)*50, 
                      shrinkage = c(0.1,0.01, 1, 10, 100),
                      n.minobsinnode = c(10,20))
  
  set.seed(567)
  FGBM <- train(response ~ ., training1, method='gbm', metric="Accuracy", trControl = numFolds, tuneGrid=grid)
  #models[[12]] <- FGBM
  
  ################ Prepare test data for validation (remove batch effect)
  val1_testing_e <- val1_testing_fix_e
  val1_testing_myeset <- val1_testing_fix_myeset
  val1_testing_pData <- val1_testing_fix_pData
  
  sample_names <- row.names(val1_testing_pData)
  val1_testing_pData$responsenstage <- paste(val1_testing_pData$Pretreatment.T.Stage ,val1_testing_pData$Pretreatment.N.Stage, sep=".")
  val1_testing_pData <- val1_testing_pData %>% 
    mutate(responsenstage = recode(responsenstage, 
                                   "1.0" = "Block1",
                                   "1.1" = "Block2",
                                   "1.2" = "Block2",
                                   "2.0" = "Block1",
                                   "2.1" = "Block2",
                                   "2.2" = "Block2",
                                   "3.0" = "Block3",
                                   "3.1" = "Block4",
                                   "3.2" = "Block4",
                                   "4.0" = "Block3",
                                   "4.1" = "Block4",
                                   "4.2" = "Block4"))
  val1_testing_pData <- as.data.frame(lapply(val1_testing_pData, function (x) if (is.factor(x)) factor(x) else x))
  #val1_testing_pData$responsetstage <- as.numeric(as.character(val1_testing_pData$responsetstage))
  phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
  
  row.names(val1_testing_pData) <- sample_names
  val1_testing_e <- exprs(val1_testing_myeset)
  
  # if(length(setdiff(row.names(edata),row.names(expression_data)))==0){
  #   combat_type <- "response+tn"
  # }else{
  #   combat_type <- "tn"
  # }
  combat_type <- "tn"
  combat <- combat_type
  if(combat=="response+tn"){
    if (!requireNamespace("sva", quietly = TRUE)){
      BiocManager::install("sva")}
    
    library(sva)
    full_em <- cbind(val1_training_e, val1_testing_e)
    full_pData <- rbind(val1_training_pData[,1:5], val1_testing_pData[,c(1:5)])
    batch <- c(rep(1, dim(val1_training_e)[2]),rep(2, dim(val1_testing_e)[2]))
    mod <- model.matrix(~as.factor(response), data=full_pData)
    set.seed(567)
    combat_edata3 <- ComBat(dat=full_em, batch=batch, mod = mod, par.prior = TRUE, ref.batch=1)
    #training_em_new <- combat_edata3[,1:dim(training_em)[2]]
    val1_testing_e <- as.matrix(combat_edata3[, (dim(val1_training_e)[2]+1):dim(full_em)[2]])
    colnames(val1_testing_e) <- row.names(val1_testing_pData)
    val1_testing_myeset <- ExpressionSet(assayData = val1_testing_e)
    phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
    ## training data her
    val1_training_e <- combat_edata3[, 1:(dim(val1_training_e)[2])]
    val1_training_myeset <- ExpressionSet(assayData = val1_training_e)
    phenoData(val1_training_myeset) <- new("AnnotatedDataFrame", data = val1_training_pData)
  }else{
    if (!requireNamespace("sva", quietly = TRUE)){
      BiocManager::install("sva")}
    
    library(sva)
    full_em <- cbind(val1_training_e, val1_testing_e)
    full_pData <- rbind(val1_training_pData[,2:4], val1_testing_pData[,c(1:3)])
    batch <- c(rep(1, dim(val1_training_e)[2]),rep(2, dim(val1_testing_e)[2]))
    #mod <- model.matrix(~as.factor(response), data=full_pData)
    set.seed(567)
    combat_edata3 <- ComBat(dat=full_em, batch=batch, par.prior = TRUE, ref.batch=1)
    #training_em_new <- combat_edata3[,1:dim(training_em)[2]]
    val1_testing_e <- as.matrix(combat_edata3[, (dim(val1_training_e)[2]+1):dim(full_em)[2]])
    colnames(val1_testing_e) <- row.names(val1_testing_pData)
    val1_testing_myeset <- ExpressionSet(assayData = val1_testing_e)
    phenoData(val1_testing_myeset) <- new("AnnotatedDataFrame", data = val1_testing_pData)
    ## training data her
    val1_training_e <- combat_edata3[, 1:(dim(val1_training_e)[2])]
    val1_training_myeset <- ExpressionSet(assayData = val1_training_e)
    phenoData(val1_training_myeset) <- new("AnnotatedDataFrame", data = val1_training_pData)
  }
  
  
  #### create 10 iterations of validation balanced set
  val1_testing_set <- list() 
  if (!requireNamespace("raster", quietly = TRUE)){
    BiocManager::install("raster")}
  set.seed(45)
  seeds <- raster::sampleInt(1000, 10, replace = FALSE)
  val_cohort_size <- "Full Cohort"
  for(iter in 1:10){
    rand <- seeds[iter]
    val1_testing_myeset_temp <- val1_testing_myeset
    val1_testing_pData_temp <- val1_testing_pData
    val1_testing_e_temp <- val1_testing_e
    testsampling <- val_cohort_size
    if(testsampling =="Balanced Cohort"){
      ## downsampling of training data
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block1")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block1")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat1 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat1 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block2")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block2")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat2 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat2 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block3")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block3")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat3 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat3 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      cat1a <- which(val1_testing_pData_temp$responsetstage=="1Block4")
      cat1b <- which(val1_testing_pData_temp$responsetstage=="0Block4")
      if(length(cat1a)< length(cat1b)){
        set.seed(rand)
        cat4 <- c(cat1a, sample(cat1b, min(length(cat1b),length(cat1a))))
      }else{
        set.seed(rand)
        cat4 <- c(cat1b, sample(cat1a, min(length(cat1b),length(cat1a))))
      }
      
      
      bal_e7 <- val1_testing_e_temp
      add_e7 <- val1_testing_e_temp[, c(cat1, cat2, cat3,cat4)]
      dim(add_e7)
      add_e7[1:5,1:5]
      dim(add_e7)
      add_e7[1:5,1:5]
      
      bal_pData <- val1_testing_pData_temp
      add_pData <- val1_testing_pData_temp[c(cat1, cat2, cat3,cat4),]
      dim(add_pData)
      add_pData[1:5, 1:5]
      
      val1_testing_pData_temp <- add_pData
      val1_testing_myeset_temp <- ExpressionSet(assayData = add_e7)
      phenoData(val1_testing_myeset_temp) <- new("AnnotatedDataFrame", data = val1_testing_pData_temp)
      val1_testing_e_temp <- exprs(val1_testing_myeset_temp)
      rm(bal_e7)
      rm(bal_pData)
      rm(add_e7)
      rm(add_pData)
      rm(cat1)
      rm(cat2)
      rm(cat3)
      rm(cat4)
      rm(cat1a)
      rm(cat1b)
    }
    ## Prepare test data for validation 1
    ## test data
    aab <- data.frame(val1_testing_pData_temp$response)
    colnames(aab) <- c("response")
    testing <- cbind(t(as.data.frame(val1_testing_e_temp)), aab)
    rm(aab)
    testing$response = factor(testing$response)
    dim(testing)
    testing$response = factor(testing$response, levels = c("responder","nonresponder"))
    val1_testing_set[[iter]] <- testing
  }
  length(intersect(row.names(val1_testing_set[[1]]), row.names(val1_testing_set[[10]])))
  
  ## check accuracies for all validation set testing combinations
  exo <-0
  acc_scores <- c()
  for(iter in 1:10){
    testing <- val1_testing_set[[iter]]
    #pro_temp <- pro+iter
    exo_temp <- exo+iter
    model <- FGBM
    bestest1 <- data.frame(significant_probes =model$coefnames)
    bestest1 <- separate_rows(bestest1, 1, sep = "`")
    bestest1 <- bestest1[bestest1$significant_probes!="",]
    testing_new <- testing[,c(bestest1$significant_probes, "response")]
    training$response = factor(training$response, levels = c("responder", "nonresponder"))
    testing_new$response = factor(testing_new$response, levels = c("responder", "nonresponder"))
    
    predictions <- model %>% predict(testing_new)
    predictions
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_new$response))){
      ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
    }else if(nlevels(as.factor(predictions))>=nlevels(as.factor(testing_new$response))){
      levels(testing_new$response) <- levels(predictions)
      ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
    }else if(nlevels(as.factor(predictions))<=nlevels(as.factor(testing_new$response))){
      levels(predictions) <- levels(testing_new$response)
      ext_accuracy <- confusionMatrix(as.factor(predictions),as.factor(testing_new$response))
    }
    acc_scores[exo_temp] <- ext_accuracy$overall[1]
  }
  
  ## Get final Model
  model <- FGBM
  testing <- val1_testing_set[[1]]
  
  bestest1 <- data.frame(significant_probes =model$coefnames)
  bestest1 <- separate_rows(bestest1, 1, sep = "`")
  bestest1 <- bestest1[bestest1$significant_probes!="",]
  testing_new <- testing[,c(bestest1$significant_probes, "response")]
  val1_final_deg_list <- deg_list[bestest1$significant_probes,]
  
  int_accuracy <- max(model$results$Accuracy, na.rm = TRUE)
  predictions <- predict(model, newdata = testing_new, type = "prob")
  row.names(predictions) <- row.names(testing_new)
  predictions_prob <- data.frame(predictions)
  predictions <- predict(model, newdata = testing_new)
  predictions
  predictions_prob$Predicted_class <- predictions
  predictions_prob$Actual_class <- testing_new$response
  
  predictions_prob1 <- predictions_prob
  thresh <- 0.50
  predictions_prob1$Predicted_class1 <- factor(ifelse(predictions_prob1$responder > thresh, "responder", "nonresponder"), levels = c("responder", "nonresponder"))
  #confusionMatrix(table(predictions_prob1$Predicted_class1, predictions_prob$Actual_class))
  
  SampleNumber <- c(1:dim(testing)[1])
  prob_plot <- predictions_prob %>% 
    ggplot(aes(x=SampleNumber, y=responder)) + 
    scale_color_manual(values=c("red", "black")) + 
    geom_point() + 
    geom_rug() + 
    ggtitle("Response scores for Testing set")+
    geom_hline(yintercept=thresh, linetype="dashed", color = "red")
  
  #rm(predictions_prob1)
  rm(thresh)
  rm(predictions_prob)
  
  ## clear
  rm(val1_testing_myeset_temp)
  rm(val1_testing_e_temp)
  rm(val1_testing_fix_pData)
  rm(full_em)
  rm(full_pData)
  rm(mod)
  rm(batch)
  rm(combat_edata3)
  rm(common)
  rm(combat)
  rm(testsampling)
  #
  return(list(prob_plot, prob_plot,prob_plot, predictions_prob1))
}

save(deg_list, public_gene_name, public_myeset, public_pData, raw_training_e, raw_training_myeset, raw_training_pData, training, bestest1,
     validate_GSE87211, predict_GSE87211, file = paste0(path_to_discovery_data,"DiscoveryBig.RData"), version = 3)


###
ex <- read.csv("Expression1.csv", header = TRUE, row.names = 1)


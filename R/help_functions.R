rm(list = ls())

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()


if(!require('magrittr'))      install.packages('magrittr')
if(!require('dplyr'))         install.packages('dplyr')
# if(!require('tidyr'))         install.packages('tidyr')
# if(!require('data.table'))    install.packages('data.table')
if(!require('reshape2'))      install.packages('reshape2')
if(!require('ggplot2'))       install.packages('ggplot2')
if(!require('tidyr'))         install.packages('tidyr')
if(!require('stringr'))       install.packages('stringr')
if(!require('optmatch'))      install.packages('optmatch')
if(!require('randomForest'))  install.packages('randomForest')
if(!require('sensitivitymw')) install.packages('sensitivitymw')

# if(!require('caret'))         install.packages('caret')


weights_fun <- function(prop, treat) {
  out <- treat + ((1 - treat) * prop / (1 - prop))
  return(out)
}


ipw.att <- function(wi, treat, yi) { 
  n <- length(treat)
  mean1 <- sum((treat == 1) * yi) / n
  mean0 <- sum((treat == 0) * wi * yi) / n
  
  out <- mean1 - mean0
  return(out)
}


prop_score_fun <- function(data, 
                           TreatCol = "Treat", 
                           OutcomeCol = "Y"
) {
  form <- as.formula(paste(TreatCol, "~." ))
  propscoremodel <- glm(form, data = data %>% 
                          dplyr::select(-OutcomeCol), family = binomial)
  propscore      <- predict(propscoremodel, type = "response") 
  return(propscore)
  
}

### IPW ATT function
ipw_att_fun <- function(data, propScores, 
                        TreatCol = "Treat", 
                        OutcomeCol = "Y") {
  data.ipw <- data.frame(data, 
                         "propscore" = propScores, 
                         "wi" = weights_fun(propScores, data[, TreatCol]))
  
  ipw.att.val <- ipw.att(wi = data.ipw$wi, treat = data.ipw[, TreatCol], yi = data.ipw[, OutcomeCol])
  
  return(ipw.att.val)
}


### T-learner

Tlearner_att_fun <- function(data, 
                             TreatCol = "Treat", 
                             OutcomeCol = "Y" ) {
  ### a model for Y(0)
  data.t0 <- data %>% 
    filter(!!sym(TreatCol) == 0) %>% 
    dplyr::select(-TreatCol)
  rf.form0 <- as.formula(paste(OutcomeCol, "~."))
  model0 <- randomForest(rf.form0, data = data.t0)
  
  
  ### Model for Y(1)
  data.t1 <- data %>% 
    filter(!!sym(TreatCol) == 1) %>% 
    dplyr::select(-TreatCol)
  rf.form1 <- as.formula(paste(OutcomeCol, "~."))
  model1 <- randomForest(rf.form1, data = data.t1)
  # model1
  
  pred1 <- predict(model1, newdata =  data)
  pred0 <- predict(model0, newdata =  data)
  
  T.att <- mean(pred1 - pred0)
  return(T.att)
}

# S-learner
Slearner_att_fun <- function(data, 
                             TreatCol = "Treat", 
                             OutcomeCol = "Y") {
  s.form <- as.formula(paste(OutcomeCol, '~.'))
  model.S <- randomForest(s.form, data = data)
  data.S1 <- data %>% 
    mutate(!!TreatCol := 1)
  data.S0 <- data %>% 
    mutate(!!TreatCol := 0)
  
  pred.S1 <- predict(model.S, newdata = data.S1)
  pred.S0 <- predict(model.S, newdata = data.S0)
  
  s.att <- mean(pred.S1 - pred.S0)
  return(s.att)
}


matching_att_fun <- function(data, 
                             TreatCol = "Treat", 
                             OutcomeCol = "Y", 
                             thresh = 1) {
  data %<>% 
    arrange(!!sym(TreatCol))
  form <- as.formula(paste(TreatCol, '~.'))
  propScore_model <- glm(form, 
                         data = data %>% 
                           dplyr::select(-OutcomeCol), 
                         family = 'binomial')
  
  if(thresh < 1) {
    data_thresh <- data[-which(propScore_model$fitted.values > thresh), ]
  } else {
    data_thresh <- data
  }
  
  propScore_model_thresh <- glm(form, 
                         data = data_thresh %>% 
                           dplyr::select(-OutcomeCol), 
                         family = 'binomial')
  distmat <- match_on(propScore_model_thresh)
  pairs <- pairmatch(distmat)
  pairs.short <- substr(pairs, start = 3, stop = 10)
  pairsnumeric <- as.numeric(pairs.short)
  tmp_tab <- cbind("index" = as.numeric(names(pairs)),
                   "pair.num" = pairsnumeric)
  index_mat <- matrix(NA, 
                      nrow = max(pairsnumeric, na.rm = T), 
                      ncol = 2)
  for(j in 1:max(pairsnumeric, na.rm = T)) {
    index_mat[j, ] <- sort(tmp_tab[which(tmp_tab[, "pair.num"] == j), "index"])
  }
  colnames(index_mat) <- c("Treat0", "Treat1")
  
  
  
  match_att <- mean(data[index_mat[, "Treat1"], OutcomeCol] -  
                      data[index_mat[, "Treat0"], OutcomeCol])
  Ttest <- t.test(data[index_mat[, "Treat1"], OutcomeCol] , 
                data[index_mat[, "Treat0"], OutcomeCol], paired = T)
  WilTest <- wilcox.test(data[index_mat[, "Treat1"], OutcomeCol] , 
                  data[index_mat[, "Treat0"], OutcomeCol], paired = T)
  data_plot <- data[c(index_mat), ] 
  return(list("match_att" = match_att, "Ttest" = Ttest, "WilTest" = WilTest, 
              "data_plot" = data_plot))

}

### Our solution: T-learner with GBM
Tgbm_learner <- function(data) {
  data <- data[sample(c(1:nrow(data)), nrow(data)), ]
  
  k <- 5 
  folds <- createFolds(c(1:nrow(data)), k)
  att_vec <- rep(NA, k)
  for (i in 1:k) {
    tmp_data <- data[-folds[[i]], ]
    data.t0 <- tmp_data %>% 
      filter(Treat == 0) %>% 
      dplyr::select(-Treat)
    rf.form0 <- as.formula(paste("Y ~."))
    model0 <- gbm(rf.form0, data = data.t0, 
                  distribution = "gaussian", 
                  shrinkage = 1, 
                  n.minobsinnode = 1, bag.fraction = 0.8, 
                  cv.folds = 10, keep.data = F, verbose = F,
                  n.trees = 200)
    
    # model0
    ### Model for Y(1)
    data.t1 <- tmp_data %>% 
      filter(Treat == 1) %>% 
      dplyr::select(-Treat)
    rf.form1 <- as.formula(paste("Y ~."))
    model1 <- gbm(rf.form1, data = data.t1, 
                  distribution = "gaussian", 
                  shrinkage = 1, 
                  n.minobsinnode = 1, bag.fraction = 0.8, 
                  cv.folds = 10, keep.data = F, verbose = F,
                  n.trees = 200)
    # model1
    # models prediction
    pred1 <- predict(model1, newdata =  tmp_data)
    pred0 <- predict(model0, newdata =  tmp_data)
    
    att_vec[i] <- mean(pred1 - pred0)
    
  } 
  return(att_vec)
}


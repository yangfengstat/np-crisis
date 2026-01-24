# 1. Library and Structure Setup
# setwd("~/Box/NP Crisis/Codes/Code/R/Reserves")
# 1.1 Library setup
library(foreign)
library(e1071)
library(imputeMissings)
library(pre)
library(nproc)
library(readxl)
library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)
#library(rfUtilities) 
# ---- Customized Functions ----

npc.signalextraction <- function(x, y, predictor.name, alpha = 0.05, delta = 0.05, split = 1, split.ratio = 0.5, n.cores = 1, randseed = 0, warning = TRUE, ...){
  set.seed(randseed)
  ind0 = which(y == 0)
  ind1 = which(y == 1)
  n0 = length(ind0)
  n1 = length(ind1)
  n0.1 = round(n0 * split.ratio)
  n1.1 = round(n1 * split.ratio)
  ind01.mat = sapply(1:split, function(i) sample(ind0, n0.1))
  ind11.mat = sapply(1:split, function(i) sample(ind1, n1.1))
  ind02.mat = sapply(1:split, function(i) setdiff(ind0, ind01.mat[, i]))
  ind12.mat = sapply(1:split, function(i) setdiff(ind1, ind11.mat[, i]))
  n0.cores = max(1, floor(n.cores/split))
  cutoffs.weights = parallel::mclapply(1:split, function(i) {
    set.seed(i + randseed)
    npc.signalextraction.split(x, y, alpha, delta, ind01.mat[, i], ind02.mat[, i], ind1, NULL, n.cores = n0.cores, warning = warning, predictor.name = predictor.name, ...)
  }, mc.cores = n.cores)
  Output <- data.frame(matrix(unlist(cutoffs.weights), nrow = 2, byrow = F), row.names = c("cutoff", "weight"))
  names(Output) <- predictor.name
  return(Output)
}

npc.signalextraction.split <- function(x, y, alpha, delta, ind01, ind02, ind11, ind12, n.cores = 1, warning = TRUE, predictor.name, ...) {
  if (missing(predictor.name) || length(predictor.name) == 0) {
    stop("predictor.name must be provided and non-empty")
  }
  train.ind = c(ind01, ind11)
  test.ind = c(ind02, ind12)
  train.x = x[train.ind, ]
  train.y = y[train.ind]
  test.x = x[test.ind, ]
  test.y = y[test.ind]
  Output <- data.frame(matrix(rep(0, NROW(predictor.name)), nrow = 2, ncol = NROW(predictor.name)), row.names = c("cutoff", "weight"))
  i <- 0
  crisis <- sum(y == 0)
  noncrisis <- sum(y == 1)
  for (predictor in predictor.name){
    i <- i+1
    obj = npc_aux.core(test.y, test.x[,predictor.name[i]], alpha = alpha, delta = delta, n.cores = n.cores, warning = warning, ...)
    Output["cutoff",i] <- obj$cutoff
    pred.label <- x[,predictor.name[i]] > Output["cutoff",i]
    missed.crisis <- if (crisis == 0) 0 else sum(y == 0 & pred.label == 1)/crisis
    false.alarm <- if (noncrisis == 0) 0 else sum(y == 1 & pred.label == 0)/noncrisis
    error <- missed.crisis + false.alarm
    
    epsilon <- 1e-9 # A small number to prevent division by zero
    weight_val <- 1 / (false.alarm + epsilon) - 1
    
    # Ensure weight is not negative
    Output["weight", i] <- max(0, weight_val)
  }
  weight_sum <- sum(Output["weight",])
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    Output["weight",] <- rep(1 / NROW(predictor.name), NROW(predictor.name))
  } else {
    Output["weight",] <- Output["weight",]/weight_sum
  }
  return(Output)
}

npc_aux.core <- function(y, score, alpha = 0.05, delta = 0.05, n.cores = 1, warning = TRUE, ...) {
  ind0 = which(y == 0)
  ind1 = which(y == 1)
  sig = TRUE
  score0 = score[ind0]
  score1 = score[ind1]
  nsmall = FALSE
  n <- length(score0)
  score.sort <- sort(score0)
  s <- 1-pbinom(0:(n-1),size=n,prob=1-alpha)
  
  valid_indices <- which(s <= delta)
  if (length(valid_indices) == 0) {
    # Handle the error: maybe return NA or a default cutoff
    if (warning) {
      base::warning("No suitable cutoff found for the given alpha/delta.")
    }
    cutoff <- NA
    loc <- NA
  } else {
    loc <- min(valid_indices)
    cutoff <- score.sort[loc]
  }
  obj <- list(scores = score0, beta.l = NULL, beta.u = NULL, alpha.l = NULL, alpha.u = NULL, ru1=NULL,rl1=NULL,ru0=NULL,rl0=NULL)
  return(list(cutoff = cutoff, loc = loc, sign = sig, alpha.l = obj$alpha.l, alpha.u = obj$alpha.u,
              beta.l = obj$beta.l, beta.u = obj$beta.u, nsmall = nsmall))
}

predict.signalextraction <- function(newx, model, predictor.name){
  cutoff <- model["cutoff",]
  weight <- model["weight",]
  flag <- data.frame(matrix(rep(0, NROW(predictor.name)), nrow = NROW(newx), ncol = NROW(predictor.name)))
  names(flag) <- predictor.name
  for (predictor in predictor.name){
    flag[,predictor] <- newx[,predictor] > cutoff[,predictor]
  }
  score <- data.matrix(flag)%*%t(data.matrix(weight))
  return(score)
}

#### Sum of Errors Calculation #####
SumErrors.label <- function(predicted, actual, weight = NULL){
  # number of crisis
  crisis <- sum(actual == 1)
  # number of non-crisis
  noncrisis <- sum(actual == 0)

  if (crisis == 0){
    missed.min <- 0
  }else{
    missed.min <- sum(predicted == 0 & actual == 1)/crisis
  }

  if (noncrisis == 0){
    false.min <- 0
  }else{
    false.min <- sum(predicted == 1 & actual == 0)/noncrisis
  }

  if (is.null(weight)){
    error.min <- false.min + missed.min
  }else{
    error.min <- weight*false.min + (1-weight)*missed.min
  }

  return(list(error = error.min, false = false.min, missed = missed.min))
}

##### Performance #####

Performance.evaluation <- function(model.trained, trainingSample.NP, testSample.NP, predictor.name){
  ##### construct temp preformance data.frame #####
  Performance.temp <- data.frame(matrix(, nrow = 1, ncol = 56))
  # Name variables in Performance
  names(Performance.temp)<-c("Error.training", "False.training", "Missed.training",
                             "Error.tuning.mean", "Error.tuning.sd",
                             "False.tuning.mean", "False.tuning.sd",
                             "Missed.tuning.mean", "Missed.tuning.sd",
                             "Error.tuning_out.mean", "Error.tuning_out.sd",
                             "False.tuning_out.mean", "False.tuning_out.sd",
                             "Missed.tuning_out.mean", "Missed.tuning_out.sd",
                             "Error.test", "False.test", "Missed.test",
                             "Error.tuning.threshold.mean", "False.tuning.threshold.mean", "Missed.tuning.threshold.mean",
                             "Error.tuning.threshold.sd", "False.tuning.threshold.sd", "Missed.tuning.threshold.sd", "Threshold.tuning.mean",
                             "Error.tuning_out.threshold.mean", "False.tuning_out.threshold.mean", "Missed.tuning_out.threshold.mean",
                             "Error.tuning_out.threshold.sd", "False.tuning_out.threshold.sd", "Missed.tuning_out.threshold.sd",
                             "Error.training.threshold", "False.training.threshold", "Missed.training.threshold", "Threshold.training",
                             "Error.test.threshold", "False.test.threshold", "Missed.test.threshold",
                             "Error.test.threshold_tuning", "False.test.threshold_tuning", "Missed.test.threshold_tuning",
                             "AUC.tuning.mean", "AUC.tuning_out.mean", "AUC.training", "AUC.test",
                             "Predictor.1st", "Predictor.2nd", "Predictor.3rd", "Predictor.4th", "Predictor.5th",
                             "Predictor.6th", "Predictor.7th", "Predictor.8th", "Predictor.9th", "Predictor.10th",
                             "Parameters")

  ##### Performance main #####
  score.training <- predict.signalextraction(trainingSample.NP, model.trained, predictor.name)
  pred.training <- score.training > 0.5

  Error.training <- SumErrors.label(pred.training, trainingSample.NP$outcome)

  Performance.temp$Error.training <- Error.training$error
  Performance.temp$False.training <- Error.training$missed
  Performance.temp$Missed.training <- Error.training$false

  ###################################### test sample performance ----
  test_size <- NROW(testSample.NP)
  if (test_size != 0){
    score.test <- predict.signalextraction(testSample.NP, model.trained, predictor.name)
    pred.test <- score.test > 0.5

    Error.test <- SumErrors.label(pred.test, testSample.NP$outcome)

    Performance.temp$Error.test <- Error.test$error
    Performance.temp$False.test <- Error.test$missed
    Performance.temp$Missed.test <- Error.test$false
  }

  return(Performance.temp)
}

rho.star <- function(pi, lambda, pi.mean, gamma, g, delta, r, sigma){
  p <- 1-delta/((1-pi.mean)*(pi.mean+delta))
  #p <- (delta * pi.mean) / ((1 - pi.mean) * (1 - delta))
  eta <- ((pi*(1-pi.mean))/(pi.mean*(1-pi))*p)^(1/sigma)
  (1 - (r-g)/(1+g)*lambda - 1/eta*(1-gamma-(1+r)/(1+g)*lambda))/(1-(1-1/eta)*delta)
}


#x: the matrix containing the predictors
#y: the response
#N: the number of slicing schemes. If missing, the slicing schemes with G=3, ..., [log(n)]+1 is used
#return: a vector containing K^G.
k.filter<-function(x,y,N=NULL,nslices=NULL,slicing.scheme=NULL,response.type="continuous",method="fused"){
  method=match.arg(arg=method,choices=c("fused","single"))
  response.type=match.arg(arg=response.type,choices=c("continuous","discrete","categorical"))

  if(!is.null(slicing.scheme)){
    if(method=="fused")warning("A slicing scheme is given. Using a single Kolmogorov filter.")
    obj<-k.filter.single(x=x,y=y,slicing.scheme=slicing.scheme)
  }else{
    if(response.type=="categorical"){
      y<-factor(y)
      obj<-k.filter.single(x=x,y=y)
    }
    if(response.type=="continuous"|response.type=="discrete"){
      n<-nrow(x)
      if(is.null(N)){N<-ceiling(log(n))-2}
      if(is.null(nslices)){nslices<-3:(N+2)}
      obj<-fused.k.filter(x=x,y=y,N=N,nslices=nslices)
    }
  }
  obj
}


k.filter.single<-function(x,y,slicing.scheme=NULL,nslices=NULL){
  if(is.factor(y)){
    y.dm<-y
  }
  if(!is.factor(y)){
    if(is.null(slicing.scheme)){
      if(is.null(nslices)){
        stop("If y is not a factor, either slicing.scheme or nslices should be specified")
      }
      slicing.scheme=quantile(quantile(y,seq(1/nslices,1-1/nslices,1/nslices)))

    }
    K<-length(slicing.scheme)+1
    slicing.scheme=c(min(y)-1,slicing.scheme,max(y)+1)
    y.dm<-cut(y,slicing.scheme,label=1:K,right=F)
  }
  y.dm<-as.numeric(y.dm)
  K<-length(unique(y.dm))
  p<-ncol(x)
  ks.stat<-matrix(0,p,K*(K-1)/2)

  nclass<-0
  for(j in 1:(K-1)){
    for(l in (j+1):K){
      nclass<-nclass+1
      for(i in 1:p){
        ks.stat[i,nclass]<-ks.test(x[y.dm==j,i],x[y.dm==l,i])$statistic}
    }
  }
  ks.stat.max0<-apply(ks.stat,1,max)
  k.rank<-rank(-ks.stat.max0,ties.method="max")
  list(k.stat=ks.stat.max0,k.rank=k.rank)
}


fused.k.filter<-function(x,y,N=NULL,nslices=NULL){
  n<-nrow(x)
  p<-ncol(x)

  if(is.null(N)){N<-ceiling(log(n))-2}
  if(is.null(nslices)){nslices<-3:(N+2)}

  ks.stat.single<-matrix(0,N,p)
  ks.stat.max<-rep(0,p)

  for(K in nslices){
    slicing.scheme<-quantile(y,seq(0,1,1/K))
    slicing.scheme[1]<-slicing.scheme[1]-1
    slicing.scheme[K+1]<-slicing.scheme[K+1]+1
    y.dm<-cut(y,slicing.scheme,labels=c(1:K),right=F)
    ks.stat<-matrix(0,p,K*(K-1)/2)

    nclass<-0
    for(j in 1:(K-1)){
      for(l in (j+1):K){
        nclass<-nclass+1
        for(i in 1:p){
          ks.stat[i,nclass]<-ks.test(x[y.dm==j,i],x[y.dm==l,i])$statistic}
      }
    }
    ks.stat.max0<-apply(ks.stat,1,max)
    ks.stat.single[K-2,]<-ks.stat.max0
    ks.stat.max<-ks.stat.max+ks.stat.max0}

  k.rank=rank(-ks.stat.max,ties.method="max")

  list(k.stat=ks.stat.max,k.stat.single=ks.stat.single,N=N,nslices=nslices,k.rank=k.rank)

}


# --- NEW: helper to fit/apply isotonic probability calibration ---
# Why: rfUtilities::probability.calibration returns calibrated probabilities but not a reusable model.
#      We learn the mapping on training, then apply it to new scores via linear interpolation.
# --- NEW: helper to fit/apply isotonic probability calibration ---
# Robust against NA/Inf, degenerate labels, or constant probabilities.
.fit_prob_calibrator <- function(y, p, regularization = TRUE, clip_eps = 1e-6) {
  if (length(y) != length(p)) stop("y and p must have the same length")
  y <- as.integer(y)
  p <- as.numeric(p)
  # basic cleaning
  p <- pmin(pmax(p, clip_eps), 1 - clip_eps)
  keep <- is.finite(p) & !is.na(y) & (y %in% c(0L, 1L))
  y <- y[keep]
  p <- p[keep]
  # guard: need both classes and at least two unique probabilities
  if (length(unique(y)) < 2 || length(unique(p)) < 2 || length(y) < 5) {
    # identity mapping
    return(list(map = function(z) pmin(pmax(as.numeric(z), 0), 1), identity = TRUE))
  }
  # try/catch around rfUtilities::probability.calibration
  p_cal <- tryCatch(
    rfUtilities::probability.calibration(y = y, p = p, regularization = regularization),
    error = function(e) {
      # fallback to identity on failure (e.g., isoreg complains about NA)
      return(NULL)
    }
  )
  if (is.null(p_cal)) {
    return(list(map = function(z) pmin(pmax(as.numeric(z), 0), 1), identity = TRUE))
  }
  # build monotone mapping p -> p_cal using sorted unique knots
  ord <- order(p)
  x <- p[ord]
  yhat <- as.numeric(p_cal[ord])
  # collapse duplicates
  knots <- aggregate(yhat ~ x, FUN = mean)
  # ensure strictly increasing x for approxfun
  ux <- unique(knots$x)
  uy <- tapply(knots$yhat, factor(match(knots$x, ux), levels = seq_along(ux)), mean)
  f <- stats::approxfun(ux, as.numeric(uy), method = "linear", rule = 2, ties = "ordered")
  list(map = f, identity = FALSE)
}

.apply_prob_calibrator <- function(calibrator, p_new) {
  p_new <- as.numeric(p_new)
  p_cal <- calibrator$map(p_new)
  # clip to [0,1]
  pmin(pmax(p_cal, 0), 1)
}

reserves_match.mse.alpha <- function(alpha.choice, current, year.cutoff_choice, test.size, Country.list, CrossSectional_Calibration, model.choice, i.model, delta = 0.1, calibrate.prob = TRUE, calibrate.regularization = TRUE){
  sigma.t <- 2
  CrossSectional_Reserves_Match <- data.frame(matrix(, nrow = length(Country.list)*length(year.cutoff_choice)*length(alpha.choice), ncol = 8))
  names(CrossSectional_Reserves_Match) <- c("Country.Name", "Year.cutoff", "Alpha", "Probability.estimate", "Reserves.observed", "Reserves.optimal", "Reserves.error", "Reserves.mse.three_years")
  country_year_alpha.id <- 0
  for (alpha.t in alpha.choice){
    for (year.cutoff in year.cutoff_choice){
      # Filter data by year and remove NA outcomes
current_working <- current[which(current$Year <= year.cutoff + test.size), ]
current_working.training <- current[which(current$Year <= year.cutoff), ]
current_working.training <- current_working.training[!is.na(current_working.training$outcome), ]

# Remove columns and rows with all NAs
      temp0 <- colSums(!is.na(current_working.training))
      current_working <- current_working[,which(temp0!=0)]
      temp0 <- rowSums(!is.na(current_working[,c(firstvar:NCOL(current_working)-1)]))
      current_working <- current_working[which(temp0!=0),]

      # Split into training and forecasting sets
      current.training <- current_working[which(current_working$Year <= year.cutoff),]


      current.forecasting <- current_working[which(current_working$Year > year.cutoff & current_working$Year <= year.cutoff+ test.size),]
      current.forecasting <- current.forecasting[!is.na(current.forecasting$outcome),]
      forecasting_row <- nrow(current.forecasting)
      current.training <- current.training[!is.na(current.training$outcome),]

      # Get predictor names excluding reserve-related variables
      predictor.name <- colnames(current.training[,c(firstvar:(NCOL(current.training)-1))])
      predictor.name <- predictor.name[which(predictor.name != "RES_BMGS" & predictor.name != "RES_BM" & predictor.name != "RES_GDP" & predictor.name != "RES_CH_GDP" & predictor.name != "RES_D_S")]
      predictor.name <- predictor.name[which(predictor.name != "DB_GDP_BIS_FC")]

      # Impute missing values using median/mode
      training.median <- imputeMissings::compute(current.training[, predictor.name], method = "median/mode")

      # Prepare training sample
      trainingSample <- current.training
      trainingSample[, predictor.name] <- imputeMissings::impute(trainingSample[, predictor.name], object  = training.median)

      # Prepare forecasting sample
      forecastingSample <- current.forecasting
      forecastingSample[, predictor.name] <- imputeMissings::impute(forecastingSample[, predictor.name], object = training.median)

      # Create NP samples with inverted outcomes
      trainingSample.NP <- trainingSample
      forecastingSample.NP <- forecastingSample
      trainingSample.NP$outcome <- 1 - trainingSample.NP$outcome
      forecastingSample.NP$outcome <- 1 - forecastingSample.NP$outcome
      if(i.model == 1){
        model.trained <- npc.signalextraction(trainingSample.NP, trainingSample.NP$outcome, predictor.name, alpha = alpha.t, delta = delta, split = 1, split.ratio = 0.5, n.cores = 1, randseed = 0)
        score.training <- data.frame(Score = 1 - as.numeric(predict.signalextraction(trainingSample.NP, model.trained, predictor.name)), Label = ifelse(1 - as.numeric(predict.signalextraction(trainingSample.NP, model.trained, predictor.name)) > 0.5, 1, 0))
        score.forecasting <- data.frame(Score = 1 - as.numeric(predict.signalextraction(forecastingSample.NP, model.trained, predictor.name)), Label = ifelse(1 - as.numeric(predict.signalextraction(forecastingSample.NP, model.trained, predictor.name)) > 0.5, 1, 0))
      } else if (i.model == 2){
        model.trained <- npc(trainingSample.NP[, predictor.name], trainingSample.NP$outcome, method = model.choice[i.model], alpha = alpha.t, delta = delta, split = 1, split.ratio = 0.5, n.cores = 1, randSeed = 0)
        score.training <- data.frame(Score = 1 - as.numeric(predict(model.trained, trainingSample.NP[, predictor.name])$pred.score), Label = 1 - as.numeric(predict(model.trained, trainingSample.NP[, predictor.name])$pred.label))
        score.forecasting <- data.frame(Score = 1 - as.numeric(predict(model.trained, forecastingSample.NP[, predictor.name])$pred.score), Label = 1 - as.numeric(predict(model.trained, forecastingSample.NP[, predictor.name])$pred.label))
      }
      
      # --- NEW: probability calibration (isotonic) ---
      # Why: improve probabilistic reliability used downstream by rho.star.
      if (calibrate.prob) {
        # Fit calibrator with training labels aligned to Score orientation (original trainingSample$outcome)
        calibrator <- .fit_prob_calibrator(y = trainingSample$outcome, p = score.training$Score, regularization = calibrate.regularization)
        # Apply to both training and forecasting probabilities
        score.training$Score   <- .apply_prob_calibrator(calibrator, score.training$Score)
        score.forecasting$Score <- .apply_prob_calibrator(calibrator, score.forecasting$Score)
      }
      # Threshold to labels for reporting paths that still need labels
      #label.training <- data.frame(Label = ifelse(score.training$Score > 0.5, 1, 0))
     # label.forecasting <- data.frame(Label = ifelse(score.forecasting$Score > 0.5, 1, 0))
      
      #posterior.training <- cbind(trainingSample[,c("Country_Code", "Country_Name", "Year", "SuddenStop_GI")], score.training)
      #posterior.training <- cbind(posterior.training, label.training)
      posterior.forecasting <- cbind(forecastingSample[,c("Country_Code", "Country_Name", "Year", "SuddenStop_GI")], score.forecasting)
    #  posterior <- rbind(posterior.training, posterior.forecasting)
      
      for (country in Country.list){
        country_year_alpha.id <- country_year_alpha.id + 1
        CrossSectional_Reserves_Match$Alpha[country_year_alpha.id] <- alpha.t
        CrossSectional_Reserves_Match$Country.Name[country_year_alpha.id] <- country
        CrossSectional_Reserves_Match$Year.cutoff[country_year_alpha.id] <- year.cutoff+1
        r <- 0.05
       if(i.model == 1){
         probability.estimate <- posterior.forecasting$Score
       }else{
         probability.estimate <- posterior.forecasting$Label
       }
        
        #probability.estimate <- posterior.forecasting$Score
        if (country == "Aggregate"){
          probability.estimate[probability.estimate == 0] <- 10^(-3)
          probability.estimate[probability.estimate == 1] <- 1-10^(-3)
          if (is_empty(posterior.forecasting$Score[which(posterior.forecasting$Score>0.001)]) | i.model == 1){
            CrossSectional_Reserves_Match$Probability.estimate[country_year_alpha.id] <- mean(probability.estimate, na.rm = TRUE)
          }else{
            CrossSectional_Reserves_Match$Probability.estimate[country_year_alpha.id] <- mean(probability.estimate, na.rm = TRUE)
          }
          CrossSectional_Reserves_Match$Reserves.observed[country_year_alpha.id] <- median(current.forecasting$RES_GDP, na.rm = TRUE)
          cur_loc <- which(CrossSectional_Calibration$Country.Name == country & CrossSectional_Calibration$Year.cutoff == year.cutoff)
          CrossSectional_Reserves_Match$Reserves.optimal[country_year_alpha.id] <- rho.star(CrossSectional_Reserves_Match$Probability.estimate[country_year_alpha.id],  
                                                                                            CrossSectional_Calibration$lambda[cur_loc], 
                                                                                            CrossSectional_Calibration$pi.mean[cur_loc], 
                                                                                            CrossSectional_Calibration$gamma[cur_loc],  
                                                                                            CrossSectional_Calibration$g[cur_loc], 
                                                                                            CrossSectional_Calibration$delta[cur_loc], 
                                                                                            r, sigma.t)*100
        }else{
          probability.estimate <- posterior.forecasting$Score[which(posterior.forecasting$Country_Name == country)]
          probability.estimate[probability.estimate == 0] <- 10^(-3)
          probability.estimate[probability.estimate == 1] <- 1-10^(-3)
          CrossSectional_Reserves_Match$Probability.estimate[country_year_alpha.id] <- ifelse(is.null(probability.estimate), NA, probability.estimate)
          CrossSectional_Reserves_Match$Reserves.observed[country_year_alpha.id] <- ifelse(is.null(current.forecasting$RES_GDP[which(current.forecasting$Country_Name == country)]), NA, current.forecasting$RES_GDP[which(current.forecasting$Country_Name == country)])
          CrossSectional_Reserves_Match$Reserves.optimal[country_year_alpha.id] <- rho.star(CrossSectional_Reserves_Match$Probability.estimate[country_year_alpha.id],
                                                                                            CrossSectional_Calibration$lambda[which(CrossSectional_Calibration$Country.Name == country & CrossSectional_Calibration$Year.cutoff == year.cutoff)], 
                                                                                            CrossSectional_Calibration$pi.mean[which(CrossSectional_Calibration$Country.Name == country & CrossSectional_Calibration$Year.cutoff == year.cutoff)], 
                                                                                            CrossSectional_Calibration$gamma[which(CrossSectional_Calibration$Country.Name == country & CrossSectional_Calibration$Year.cutoff == year.cutoff)],  
                                                                                            CrossSectional_Calibration$g[which(CrossSectional_Calibration$Country.Name == country & CrossSectional_Calibration$Year.cutoff == year.cutoff)], 
                                                                                            CrossSectional_Calibration$delta[which(CrossSectional_Calibration$Country.Name == country & CrossSectional_Calibration$Year.cutoff == year.cutoff)], 
                                                                                            r, sigma.t)*100
        }
        CrossSectional_Reserves_Match$Reserves.error[country_year_alpha.id] <-  CrossSectional_Reserves_Match$Reserves.optimal[country_year_alpha.id] -  CrossSectional_Reserves_Match$Reserves.observed[country_year_alpha.id]
        CrossSectional_Reserves_Match$Reserves.error.squared[country_year_alpha.id] <-  (CrossSectional_Reserves_Match$Reserves.optimal[country_year_alpha.id] -  CrossSectional_Reserves_Match$Reserves.observed[country_year_alpha.id])^2
        #print(paste(model.choice[i.model], "for", as.character(year.cutoff), "with alpha", as.character(alpha.t), "is trained successfully.", sep=" "))
      }
    }
    print(paste(model.choice[i.model], "with alpha0 =", as.character(alpha.t), "is trained successfully.", sep=" "))
    for (country in Country.list){
      for (year.cutoff in year.cutoff_choice){
        cur_ind <- which(CrossSectional_Reserves_Match$Alpha == alpha.t & 
                           CrossSectional_Reserves_Match$Country.Name == country & 
                           CrossSectional_Reserves_Match$Year.cutoff == year.cutoff)
        cur_ind_plus_minus_1 <- which(CrossSectional_Reserves_Match$Alpha == alpha.t & 
                                                   CrossSectional_Reserves_Match$Country.Name == country & 
                                                   CrossSectional_Reserves_Match$Year.cutoff >= year.cutoff-1 & 
                                                   CrossSectional_Reserves_Match$Year.cutoff <= year.cutoff+1)
        CrossSectional_Reserves_Match$Reserves.mse.three_years[cur_ind] <- 
          mean(CrossSectional_Reserves_Match$Reserves.error[cur_ind_plus_minus_1]^2, na.rm = TRUE)
      }
    }
  }
  return(CrossSectional_Reserves_Match)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#   Functions IOTA5 analyses    #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### Packages
library(car)
library(metafor)
library(boot)
library(auRoc)
library(data.table)
library(SDMTools)
library(doBy)
library(matrixStats)
library(Metatron)
library(logistf)
library(plotrix)


###########################
### 1. Function for AUC ###
###########################

## Without multiple imputation
AUC.IOTA <- function(pred, outcome, center, data, method.MA = "REML", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  AUCcenter <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCcenter <- data.frame(AUCcenter)
  AUCcenter$Center <- centers
  
  # AUC per center
  for(i in seq_along(centers)){
    AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "pepe")
    AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i],])
    AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 1 & Df$center == centers[i],])/nrow(Df[Df$center == centers[i],])*100)
    
    ## Additional part for AUCs of 1
    if(AUCcenter[i, 4] == 1){
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "newcombe") # Newcombe ipv pepe
    } else{
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "pepe")
    }
    
    if(AUCcenter$AUC[i] != 1){
      AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
      AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    } else{
      AUCcenter$logit.AUC[i] <- logit(0.999)
      AUCcenter$logit.se[i]  <- (logit(0.999) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    }
  }
  
  
  AUCcenter$logit.AUC <- logit(AUCcenter$AUC)
  AUCcenter$logit.se  <- (logit(AUCcenter$AUC) - logit(AUCcenter$LL))/1.96
  AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df)
  AUCoverall$Prevalence <- round(nrow(Df[Df$y == 1,])/nrow(Df)*100)
  
  # Meta-analyse voor overall estimate
  fit.RE = rma.uni(AUCcenter$logit.AUC, sei = AUCcenter$logit.se, method = method.MA)
  PI = predict(fit.RE, transf = transf.ilogit)
  
  AUCoverall$AUC[1] <- inv.logit(coef(fit.RE))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- PI$pred
  AUCoverall$LL[2] <- PI$cr.lb
  AUCoverall$UL[2] <- PI$cr.ub
  AUCoverall
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcenter[, 1:6], NAforest, AUCoverall)
  
  # Layout forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUC (95% CI)', 'N')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcenter, IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}


## With multiple imputation
AUCimp.IOTA <- function(pred, outcome, center, imp, data, method.MA = "REML", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  AUCimp <- list()
  for(i in 1:NRimp){
    AUCimp[[i]] <- list()
  }
  
  PrevalenceOverall <- matrix(nrow = NRimp, ncol = 1)
  
  # AUC per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    AUCcenter <- matrix(ncol = 6, nrow = length(centers))
    colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
    AUCcenter <- data.frame(AUCcenter)
    AUCcenter$Center <- centers
    
    for(i in seq_along(centers)){
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 0 & Df$imp == j], method = "pepe")
      AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i] & Df$imp == j,])
      AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 1 & Df$center == centers[i] & Df$imp == j,])/nrow(Df[Df$center == centers[i] & Df$imp == j,])*100)
      
      ## Additional part for AUCs of 1
      if(AUCcenter[i, 4] == 1){
        AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 0 & Df$imp == j], method = "newcombe") # Newcombe ipv pepe
      } else{
        AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 0 & Df$imp == j], method = "pepe")
      }
      
      if(AUCcenter$AUC[i] != 1){
        AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
        AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
        AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
      } else{
        AUCcenter$logit.AUC[i] <- logit(0.999)
        AUCcenter$logit.se[i]  <- (logit(0.999) - logit(AUCcenter$LL[i]))/1.96
        AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
      }
      
    }
    AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
    
    AUCimp[[j]] <- AUCcenter
    
    PrevalenceOverall[j] <- round(nrow(Df[Df$y == 1 & Df$imp == j,])/nrow(Df[Df$imp == j,])*100)
  }
  
  AUCimpLong <- rbindlist(AUCimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  AUCcombined <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcombined) <- c('Center', 'SampleSize', 'Prevalence', 'logit.AUC', 'logit.LL', 'logit.UL')
  AUCcombined <- data.frame(AUCcombined)
  AUCcombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    AUCcombined$SampleSize[i] <- unique(AUCimpLong$SampleSize[AUCimpLong$Center == centers[i]])
    AUCcombined$Prevalence[i] <- round(mean(AUCimpLong$Prevalence[AUCimpLong$Center == centers[i]]))
    AUCcombined[i, 4] <- mean(AUCimpLong$logit.AUC[AUCimpLong$Center == centers[i]])
    WithinVar <- mean(AUCimpLong$logit.var[AUCimpLong$Center == centers[i]])
    BetweenVar <- var(AUCimpLong$logit.AUC[AUCimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    AUCcombined$PooledSE[i] <- sqrt(PooledVar)
    AUCcombined$logit.LL[i] <- AUCcombined$logit.AUC[i] - 1.96*AUCcombined$PooledSE[i]
    AUCcombined$logit.UL[i] <- AUCcombined$logit.AUC[i] + 1.96*AUCcombined$PooledSE[i]
  }
  
  AUCcombined$AUC <- inv.logit(AUCcombined$logit.AUC)
  AUCcombined$LL <- inv.logit(AUCcombined$logit.LL)
  AUCcombined$UL <- inv.logit(AUCcombined$logit.UL)
  AUCcombined <- AUCcombined[order(AUCcombined$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df[Df$imp == 1,])
  AUCoverall$Prevalence <- round(mean(PrevalenceOverall))
  
  # Meta-analyse voor overall estimate
  fit.RE = rma.uni(AUCcombined$logit.AUC, sei = AUCcombined$PooledSE, method = method.MA)
  PI = predict(fit.RE, transf = transf.ilogit)
  
  AUCoverall$AUC[1] <- inv.logit(coef(fit.RE))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- PI$pred
  AUCoverall$LL[2] <- PI$cr.lb
  AUCoverall$UL[2] <- PI$cr.ub
  AUCoverall
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcombined[, c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')], NAforest, NAforest, AUCoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUC (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcombined[, c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')], IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}

## Pooled estimate
AUC.imp <- function(pred, outcome, imp, data){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, imp = imp, stringsAsFactors = F)
  
  NRimp <- length(unique(Df$imp))
  
  AUCimp <- matrix(ncol = 3, nrow = NRimp)
  colnames(AUCimp) <- c('AUC', 'LL', 'UL')
  AUCimp <- data.frame(AUCimp)
  
  # AUC per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    AUC <- auc.nonpara.mw(Df$p[Df$y == 1 & Df$imp == j], Df$p[Df$y == 0 & Df$imp == j], method = "pepe")
    
    if(AUC[1] < 0.50){
      AUC2 <- auc.nonpara.mw(Df$p[Df$y == 0 & Df$imp == j], Df$p[Df$y == 1 & Df$imp == j], method = "pepe")
      AUCimp[j, 1:3] <- AUC2 # Zet output in tabel
    } 
    
    ## Additional part for AUCs of 1
    else{
      if(AUC[1] == 1){
        AUC3 <- auc.nonpara.mw(Df$p[Df$y == 1 & Df$imp == j], Df$p[Df$y == 0 & Df$imp == j], method = "newcombe") # Newcombe ipv pepe
        AUCimp[j, 1:3] <- AUC3 # Zet output in tabel
      } else{
        AUCimp[j, 1:3] <- AUC # Zet output in tabel
      }
    }
    
    if(AUCimp$AUC[j] != 1){
      AUCimp$logit.AUC[j] <- logit(AUCimp$AUC[j])
      AUCimp$logit.se[j]  <- (logit(AUCimp$AUC[j]) - logit(AUCimp$LL[j]))/1.96
      AUCimp$logit.var[j] <- AUCimp$logit.se[j]^2
    } else{
      AUCimp$logit.AUC[j] <- logit(0.999)
      AUCimp$logit.se[j]  <- (logit(0.999) - logit(AUCimp$LL[j]))/1.96
      AUCimp$logit.var[j] <- AUCimp$logit.se[j]^2
    }
    AUCimp$Mal[j] <- nrow(Df[Df$y == 1 & Df$imp == j,])
  }
  
  # Combine results with Rubin's rule
  AUCcombined <- matrix(ncol = 3, nrow = 1)
  colnames(AUCcombined) <- c('logit.AUC', 'logit.LL', 'logit.UL')
  AUCcombined <- data.frame(AUCcombined)
  
  AUCcombined$logit.AUC <- mean(AUCimp$logit.AUC)
  WithinVar <- mean(AUCimp$logit.var)
  BetweenVar <- var(AUCimp$logit.AUC)
  PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
  AUCcombined$PooledSE <- sqrt(PooledVar)
  AUCcombined$logit.LL <- AUCcombined$logit.AUC - 1.96*AUCcombined$PooledSE
  AUCcombined$logit.UL <- AUCcombined$logit.AUC + 1.96*AUCcombined$PooledSE
  
  # Transform back to original scale
  AUCcombined$AUC <- inv.logit(AUCcombined$logit.AUC)
  AUCcombined$LL <- inv.logit(AUCcombined$logit.LL)
  AUCcombined$UL <- inv.logit(AUCcombined$logit.UL)
  
  return(structure(list(Performance = AUCcombined, data = Df, imp = AUCimp)))
  
}


###################################################
### 2. Function for sensitivity and specificity ###
###################################################


## With multiple imputation: fixed threshold
SS.imp <- function(pred, outcome, threshold, center, imp, data, method.MA = "REML"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  threshold = threshold
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  Sensimp <- list()
  for(i in 1:NRimp){
    Sensimp[[i]] <- list()
  }
  
  
  # Sensitivity and specificity per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    Senscenter <- matrix(ncol = 7, nrow = length(centers))
    colnames(Senscenter) <- c("Center", "TP", "TN", "FP", "FN", "Sensitivity", "Specificity")
    Senscenter <- data.frame(Senscenter)
    Senscenter$Center <- centers
    
    for(i in seq_along(centers)){
      predicted_values <- ifelse(Df$p[Df$center == centers[i] & Df$imp == j] >= threshold, 1, 0)
      actual_values <- Df$y[Df$center == centers[i] & Df$imp == j]
      conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
      Senscenter$TN[i]     <- conf_matrix[1,1]
      Senscenter$TP[i]     <- conf_matrix[2,2]
      Senscenter$FN[i]     <- conf_matrix[1,2]
      Senscenter$FP[i]     <- conf_matrix[2,1]
      
      if(any(Senscenter[i, 2:5] == 0)){
        Senscenter$TN[i]     <- conf_matrix[1,1] + 0.1
        Senscenter$TP[i]     <- conf_matrix[2,2] + 0.1
        Senscenter$FN[i]     <- conf_matrix[1,2] + 0.1
        Senscenter$FP[i]     <- conf_matrix[2,1] + 0.1
      }
      
    }
    
    Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
    Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
    Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
    Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
    
    # Logit transformation
    Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
    Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
    Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
    
    Senscenter$logit.spec     <- logit(Senscenter$Specificity)
    Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
    Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
    
    Sensimp[[j]] <- Senscenter
    
  }
  
  SensimpLong <- rbindlist(Sensimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  Senscombined <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscombined) <- c('Center', 'logit.sens', 'logit.LL.sens', 'logit.UL.sens', 'logit.spec', 'logit.LL.spec', 'logit.UL.spec')
  Senscombined <- data.frame(Senscombined)
  Senscombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    
    Senscombined$TP[i] <- mean(SensimpLong$TP[SensimpLong$Center == centers[i]])
    Senscombined$TN[i] <- mean(SensimpLong$TN[SensimpLong$Center == centers[i]])
    Senscombined$FP[i] <- mean(SensimpLong$FP[SensimpLong$Center == centers[i]])
    Senscombined$FN[i] <- mean(SensimpLong$FN[SensimpLong$Center == centers[i]])
    
    # Sensitivity
    Senscombined$logit.sens[i] <- mean(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    WithinVar  <- mean(SensimpLong$logit.var.sens[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    PooledVar  <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.sens[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.sens[i] <- Senscombined$logit.sens[i] - 1.96*Senscombined$PooledSE.sens[i]
    Senscombined$logit.UL.sens[i] <- Senscombined$logit.sens[i] + 1.96*Senscombined$PooledSE.sens[i]
    
    # Specificity
    Senscombined$logit.spec[i] <- mean(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    WithinVar <- mean(SensimpLong$logit.var.spec[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.spec[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.spec[i] <- Senscombined$logit.spec[i] - 1.96*Senscombined$PooledSE.spec[i]
    Senscombined$logit.UL.spec[i] <- Senscombined$logit.spec[i] + 1.96*Senscombined$PooledSE.spec[i]
  }
  
  # Transform back to original scale
  Senscombined$Sensitivity <- inv.logit(Senscombined$logit.sens)
  Senscombined$LL.sens     <- inv.logit(Senscombined$logit.LL.sens)
  Senscombined$UL.sens     <- inv.logit(Senscombined$logit.UL.sens)
  
  Senscombined$Specificity <- inv.logit(Senscombined$logit.spec)
  Senscombined$LL.spec     <- inv.logit(Senscombined$logit.LL.spec)
  Senscombined$UL.spec     <- inv.logit(Senscombined$logit.UL.spec)
  
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  
  ## Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscombined[, c("Center", "logit.sens", "logit.spec", "PooledSE.sens", "PooledSE.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec-1 | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscombined, IncludedCenter = centers, data = Df)))
}


## With multiple imputation: fixed sensitivity
FixedSens.imp <- function(pred, outcome, Sensitivity, center, imp, data, method.MA = "REML"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  Sensitivity = Sensitivity
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  Sensimp <- list()
  for(i in 1:NRimp){
    Sensimp[[i]] <- list()
  }
  
  
  # Sensitivity and specificity per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    Senscenter <- matrix(ncol = 3, nrow = length(centers))
    colnames(Senscenter) <- c("Center", "Sensitivity", "Specificity")
    Senscenter <- data.frame(Senscenter)
    Senscenter$Center <- centers
    Senscenter$Sensitivity <- Sensitivity
    
    for(i in seq_along(centers)){
      threshold <- quantile(Df$p[Df$center == centers[i] & Df$imp == j & Df$y == 1], probs = 1 - Sensitivity)
      
      predicted_values <- ifelse(Df$p[Df$center == centers[i] & Df$imp == j] >= threshold, 1, 0)
      actual_values <- Df$y[Df$center == centers[i] & Df$imp == j]
      conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
      Senscenter$TN[i]     <- conf_matrix[1,1]
      Senscenter$TP[i]     <- conf_matrix[2,2]
      Senscenter$FN[i]     <- conf_matrix[1,2]
      Senscenter$FP[i]     <- conf_matrix[2,1]
      
      if(any(Senscenter[i, 4:7] == 0)){
        Senscenter$TN[i]     <- conf_matrix[1,1] + 0.1
        Senscenter$TP[i]     <- conf_matrix[2,2] + 0.1
        Senscenter$FN[i]     <- conf_matrix[1,2] + 0.1
        Senscenter$FP[i]     <- conf_matrix[2,1] + 0.1
      }
      
    }
    
    Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
    Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
    Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
    Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
    
    # Logit transformation
    Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
    Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
    Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
    
    Senscenter$logit.spec     <- logit(Senscenter$Specificity)
    Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
    Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
    
    Sensimp[[j]] <- Senscenter
    
  }
  
  SensimpLong <- rbindlist(Sensimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  Senscombined <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscombined) <- c('Center', 'logit.sens', 'logit.LL.sens', 'logit.UL.sens', 'logit.spec', 'logit.LL.spec', 'logit.UL.spec')
  Senscombined <- data.frame(Senscombined)
  Senscombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    
    Senscombined$TP[i] <- mean(SensimpLong$TP[SensimpLong$Center == centers[i]])
    Senscombined$TN[i] <- mean(SensimpLong$TN[SensimpLong$Center == centers[i]])
    Senscombined$FP[i] <- mean(SensimpLong$FP[SensimpLong$Center == centers[i]])
    Senscombined$FN[i] <- mean(SensimpLong$FN[SensimpLong$Center == centers[i]])
    
    # Sensitivity
    Senscombined$logit.sens[i] <- mean(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    WithinVar  <- mean(SensimpLong$logit.var.sens[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    PooledVar  <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.sens[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.sens[i] <- Senscombined$logit.sens[i] - 1.96*Senscombined$PooledSE.sens[i]
    Senscombined$logit.UL.sens[i] <- Senscombined$logit.sens[i] + 1.96*Senscombined$PooledSE.sens[i]
    
    # Specificity
    Senscombined$logit.spec[i] <- mean(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    WithinVar <- mean(SensimpLong$logit.var.spec[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.spec[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.spec[i] <- Senscombined$logit.spec[i] - 1.96*Senscombined$PooledSE.spec[i]
    Senscombined$logit.UL.spec[i] <- Senscombined$logit.spec[i] + 1.96*Senscombined$PooledSE.spec[i]
  }
  
  # Transform back to original scale
  Senscombined$Sensitivity <- inv.logit(Senscombined$logit.sens)
  Senscombined$LL.sens     <- inv.logit(Senscombined$logit.LL.sens)
  Senscombined$UL.sens     <- inv.logit(Senscombined$logit.UL.sens)
  
  Senscombined$Specificity <- inv.logit(Senscombined$logit.spec)
  Senscombined$LL.spec     <- inv.logit(Senscombined$logit.LL.spec)
  Senscombined$UL.spec     <- inv.logit(Senscombined$logit.UL.spec)
  
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  

  ## Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscombined[, c("Center", "logit.sens", "logit.spec", "PooledSE.sens", "PooledSE.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscombined, IncludedCenter = centers, data = Df)))
}


## With multiple imputation: fixed specificity
FixedSpec.imp <- function(pred, outcome, Specificity, center, imp, data, method.MA = "REML"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  Specificity = Specificity
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  Sensimp <- list()
  for(i in 1:NRimp){
    Sensimp[[i]] <- list()
  }
  
  
  # Sensitivity and specificity per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    Senscenter <- matrix(ncol = 3, nrow = length(centers))
    colnames(Senscenter) <- c("Center", "Sensitivity", "Specificity")
    Senscenter <- data.frame(Senscenter)
    Senscenter$Center <- centers
    Senscenter$Specificity <- Specificity
    
    for(i in seq_along(centers)){
      threshold <- quantile(Df$p[Df$center == centers[i] & Df$imp == j & Df$y == 0], probs = Specificity)
      
      predicted_values <- ifelse(Df$p[Df$center == centers[i] & Df$imp == j] >= threshold, 1, 0)
      actual_values <- Df$y[Df$center == centers[i] & Df$imp == j]
      conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
      Senscenter$TN[i]     <- conf_matrix[1,1]
      Senscenter$TP[i]     <- conf_matrix[2,2]
      Senscenter$FN[i]     <- conf_matrix[1,2]
      Senscenter$FP[i]     <- conf_matrix[2,1]
      
      if(any(Senscenter[i, 4:7] == 0)){
        Senscenter$TN[i]     <- conf_matrix[1,1] + 0.1
        Senscenter$TP[i]     <- conf_matrix[2,2] + 0.1
        Senscenter$FN[i]     <- conf_matrix[1,2] + 0.1
        Senscenter$FP[i]     <- conf_matrix[2,1] + 0.1
      }
      
    }
    
    Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
    Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
    Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
    Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
    
    # Logit transformation
    Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
    Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
    Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
    
    Senscenter$logit.spec     <- logit(Senscenter$Specificity)
    Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
    Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
    
    Sensimp[[j]] <- Senscenter
    
  }
  
  SensimpLong <- rbindlist(Sensimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  Senscombined <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscombined) <- c('Center', 'logit.sens', 'logit.LL.sens', 'logit.UL.sens', 'logit.spec', 'logit.LL.spec', 'logit.UL.spec')
  Senscombined <- data.frame(Senscombined)
  Senscombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    
    Senscombined$TP[i] <- mean(SensimpLong$TP[SensimpLong$Center == centers[i]])
    Senscombined$TN[i] <- mean(SensimpLong$TN[SensimpLong$Center == centers[i]])
    Senscombined$FP[i] <- mean(SensimpLong$FP[SensimpLong$Center == centers[i]])
    Senscombined$FN[i] <- mean(SensimpLong$FN[SensimpLong$Center == centers[i]])
    
    # Sensitivity
    Senscombined$logit.sens[i] <- mean(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    WithinVar  <- mean(SensimpLong$logit.var.sens[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    PooledVar  <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.sens[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.sens[i] <- Senscombined$logit.sens[i] - 1.96*Senscombined$PooledSE.sens[i]
    Senscombined$logit.UL.sens[i] <- Senscombined$logit.sens[i] + 1.96*Senscombined$PooledSE.sens[i]
    
    # Specificity
    Senscombined$logit.spec[i] <- mean(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    WithinVar <- mean(SensimpLong$logit.var.spec[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.spec[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.spec[i] <- Senscombined$logit.spec[i] - 1.96*Senscombined$PooledSE.spec[i]
    Senscombined$logit.UL.spec[i] <- Senscombined$logit.spec[i] + 1.96*Senscombined$PooledSE.spec[i]
  }
  
  # Transform back to original scale
  Senscombined$Sensitivity <- inv.logit(Senscombined$logit.sens)
  Senscombined$LL.sens     <- inv.logit(Senscombined$logit.LL.sens)
  Senscombined$UL.sens     <- inv.logit(Senscombined$logit.UL.sens)
  
  Senscombined$Specificity <- inv.logit(Senscombined$logit.spec)
  Senscombined$LL.spec     <- inv.logit(Senscombined$logit.LL.spec)
  Senscombined$UL.spec     <- inv.logit(Senscombined$logit.UL.spec)
  
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"

  
  ## Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscombined[, c("Center", "logit.sens", "logit.spec", "PooledSE.sens", "PooledSE.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscombined[, c(1, 10:15)], IncludedCenter = centers, data = Df)))
}


###################################
### 3. Function for calibration ###
###################################

## RMI: Pooled data
CalPool.RMI <- function(pred, outcome, imp, data,
                        title = "Calibration curve", xlab = "Risk of Malignancy Index (RMI)", 
                        ylab = "Observed proportion of malignancy"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  LP = log(as.numeric(eval(arguments$pred, data)) + 1)
  outcome <- eval(arguments$outcome, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, LP = LP, imp = imp, stringsAsFactors = F)
  
  NrImp = length(unique(imp))
  
  # Matrices
  FE <- matrix(ncol = 3, nrow = NrImp)
  colnames(FE) <- c("Imp", "Intercept", "LP")
  FE <- data.frame(FE)
  FE$Imp <- seq(1, NrImp, by = 1)
  
  RE <- list()
  for(i in 1:NrImp){
    RE[[i]] <- list()
  }
  

  calplotdata.avg <- as.data.frame(matrix(NA, NrImp, 2*NrImp + 1))
  colnames(calplotdata.avg) <-  c("LP", paste("obs", c(1:NrImp)), paste("var", c(1:NrImp)))
  calplotdata.avg$LP <-  seq(min(Df$LP), log(1001), length.out = NrImp)
  
  
  for(i in 1:NrImp){
    cat("Imputation", i, " of", NrImp, "\n")
    cal.model <- glm(y ~ 1 + LP, data = Df[Df$imp == i, ], family = "binomial")
    FE$Intercept[i] <- coef(cal.model)[1]
    FE$LP[i] <- coef(cal.model)[2]

    # Predict probability
    calplotdata.avg[, i + 1] <- predict(cal.model, newdata = calplotdata.avg, type = "response", allow.new.levels = T)
    
    #get the variance for the confidence interval
    mm <- model.matrix( ~ calplotdata.avg$LP, calplotdata.avg)
    calplotdata.avg[, i + NrImp + 1] <-
      diag(mm %*% tcrossprod(vcov(cal.model), mm))
  }  
  
  calplotdata.avg$p_pred <- exp(calplotdata.avg$LP) - 1
  calplotdata.avg$p_obs <- rowMeans(calplotdata.avg[ ,2:NrImp + 1], na.rm=T)
  calplotdata.avg$lp_obs <- rowMeans(calplotdata.avg[ ,2:NrImp + 1], na.rm=T)
  
  #combine between and within imputatation variance
  Bvar <- NrImp + 2
  Evar <- 2 * NrImp + 1
  Eobs <- NrImp + 1
  calplotdata.avg$var.w<-rowMeans(calplotdata.avg[ ,Bvar:Evar], na.rm=T) #gemiddelde variance binnen imputaties
  calplotdata.avg$var.b<-rowVars(as.matrix(calplotdata.avg[ ,2:Eobs]), na.rm=T) #between-imputation variance
  calplotdata.avg$var.t<-with(calplotdata.avg, var.w+(1+1/NrImp)*var.b)
  calplotdata.avg$se.t<-with(calplotdata.avg, sqrt(var.t))
  calplotdata.avg$df<-with(calplotdata.avg, (1+(NrImp/(1+NrImp))*var.w/var.b)^2*(NrImp-1))
  calplotdata.avg$t<-with(calplotdata.avg, qt(.975, df))
  calplotdata.avg$lcl<-with(calplotdata.avg, lp_obs-t*se.t)
  calplotdata.avg$ucl<-with(calplotdata.avg, lp_obs+t*se.t)
  calplotdata.avg$p_plo <- plogis(calplotdata.avg$lcl)
  calplotdata.avg$p_phi <- plogis(calplotdata.avg$ucl)
  
  par(pty = 's')
  
  plot(y = calplotdata.avg$p_obs,
       x = calplotdata.avg$p_pred,
       type = "l",
       las = 1,
       xlab = xlab,
       ylab = ylab,
       main = title,
       col = "red",
       ylim = c(0, 1),
       lwd = 1,
       cex.lab = 1, # cex.lab = 1.5
       cex.axis = 1 # cex.axis = 1.5
  )
  lines(calplotdata.avg$p_pred,
        calplotdata.avg$p_obs,
        col = "black",
        lwd = 2)
  
  return(structure(list(Predictions = calplotdata.avg, data = Df)))
}

## RMI: Imputed data
RE.ValProbImp.RMI <- function (p, y, center, imputation.id, patientid, data, LogCal = T, 
                               flexible = F, CL = c("none", "CI", "PI"), CalibrLines = c("overall", 
                                                                                         "centers", "both"), dostats = T, statloc = c(0, 0.85), 
                               legendloc = c(500, 0.15), roundstats = 3, cex = 0.75, cex.leg = 0.75, 
                               ncol.leg = 1, lty.overall = 1, lwd.overall = 2, col.overall = "black", 
                               RMprompt = F, RmNonConv = F, title = "", xlab = "Risk of Malignancy Index (RMI)", 
                               ylab = "Observed proportion of malignancy", xlim = c(0, 1000), ylim = c(0, 
                                                                                                       1), d0lab = "0", d1lab = "1", cex.d01 = 0.7, dist.label = 0.04, 
                               line.bins = -0.05, dist.label2 = 0.03, las = 1, length.seg = 1, 
                               y.intersp = 1, #lty.ideal = 1, col.ideal = "red", lwd.ideal = 1, 
                               lty.centers = NULL, lwd.centers = NULL, col.centers = NULL, 
                               Parallel = F, NrCores = detectCores() - 1, alpha = 0.05, 
                               EpsGrad = 0.001, fNonConv = c("warning", "stop"), Controlglmer = glmerControl(optimizer = "bobyqa"), 
                               AUCmeth = c("RE.auc.imp", "RE.auc.imp2"), ...)
{
  Argz = as.list(match.call())[-1]
  p = eval(Argz$p, data)
  y = eval(Argz$y, data)
  center = eval(Argz$center, data)
  id = eval(Argz$patientid, data)
  LP = log(p + 1)
  CL = match.arg(CL)
  CalibrLines = match.arg(CalibrLines)
  AUCmeth = match.arg(AUCmeth)
  if (is.factor(center)) 
    center = as.character(center)
  if (length(unique(center)) == 1) 
    stop("There is only one center, hence a random effects model should not be used here.")
  if (length(unique(center)) < 5) 
    warning("There are less than 5 centers, consider using a different method.", 
            immediate. = T)
  if (LogCal & flexible) 
    stop("LogCal and flexible cannot both be true.")
  if (missing(imputation.id)) 
    stop("Specify the variable indicating the number of the imputation!")
  imp = eval(Argz$imputation.id, data)
  if (!is.numeric(imp)) 
    imp = as.numeric(imp)
  NrImp = length(unique(imp))
  Df = data.frame(y = y, p = p, LP = log(p + 1), center = center, 
                  imp = imp, id = id)
  Df = Df[with(Df, order(center, p)), ]
  RmCenter = dlply(Df[Df$imp == 1, ], .(center), function(x) if (sum(x$y == 
                                                                     1) < 10 | sum(x$y == 0) < 10) 
    unique(x$center)
    else NULL)
  RmCenter = REMA:::ListToVector(RmCenter)
  if (length(RmCenter) != 0) {
    if (RMprompt) {
      Cmmnd = readline(paste0("The center(s) ", paste0(RmCenter, 
                                                       collapse = ", "), " have less than 10 (non-)events and these will be", 
                              " removed. Do you want to continue? (y/n)   "))
      if (Cmmnd == "n") 
        stop("Function was stopped by the user.")
    }
    Df = Df[!(Df$center %in% RmCenter), ]
  }
  IncludedCenters = unique(Df$center)
  centers = IncludedCenters
  nCenters = length(IncludedCenters)
  if (!Parallel & length(unique(Df$imp)) > 50) {
    cat("\n\nParallel computing can be used to get the results faster.\n\n")
  }
  else if (Parallel & length(unique(Df$imp)) < 10) {
    cat(paste("\n\nParallel computing is only recommended when there is a large number of imputed datasets.", 
              " When used with a small number of imputed datasets, this may be slower.\n\n"))
  }
  if (CalibrLines != "overall") {
    if (all(sapply(list(lty.centers, lwd.centers, col.centers), 
                   is.null))) {
      lty.centers = rep(1, nCenters)
      lwd.centers = rep(1, nCenters)
      col.centers = seq_len(nCenters)
    }
    else if (sum(!sapply(list(lty.centers, lwd.centers, col.centers), 
                         is.null)) != 3) {
      stop("If you specify one of the arguments, please also specify the others.")
    }
    else {
      if (any(sapply(list(lty.centers, lwd.centers, col.centers), 
                     function(x) length(x) < nCenters))) 
        stop("The vector length of lty.centers, lwd.centers and col.centers is less than the number of centers.")
      FixL = function(x) x[1:nCenters]
      lty.centers = FixL(lty.centers)
      col.centers = FixL(col.centers)
      lwd.centers = FixL(lwd.centers)
    }
  }
  RubinRes <- function(x) {
    do.call("rbind", lapply(x, function(x) {
      tmp = REMA:::rubinrules(Est, Var, x, alpha = alpha)[c("Mean estimate", 
                                                            "Total variance")]
      tmp = do.call("cbind.data.frame", tmp)
      names(tmp) = c("Est", "Var")
      tmp
    }))
  }
  RubinCL <- function(x) {
    do.call("rbind", lapply(x, function(x) {
      tmp = REMA:::rubinrules(Est, Var, x, alpha = alpha, expit = T)$Expit
      names(tmp) = c("Est", "LCL", "UCL")
      tmp
    }))
  }
  NewX = data.frame(LP = seq(min(Df$LP), max(Df$LP), length = 500))
  PerC = data.frame(LP = rep(NewX$LP, length(IncludedCenters)), 
                    center = sort(rep(IncludedCenters, 500)))
  ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
    tmp = Df[Df$imp == i, ]
    LogMM = glmer(y ~ LP + (LP | center), data = tmp, family = binomial, 
                  control = Controlglmer)
    X = if (CL != "PI") {
      as.matrix(cbind(1, NewX$LP))
    } else {
      model.matrix(LogMM)[order(tmp$p), ]
    }
    FE = fixef(LogMM)
    VarFE = diag(as.matrix(vcov.merMod(LogMM)))
    Est = X %*% FE
    Var = diag(X %*% tcrossprod(as.matrix(vcov.merMod(LogMM)), 
                                X))
    Results = list(B = list(Interc = cbind.data.frame(Est = FE[1], 
                                                      Var = VarFE[1]), Slope = cbind.data.frame(Est = FE[2], 
                                                                                                Var = VarFE[2])), LP = cbind.data.frame(Est = Est, 
                                                                                                                                        Var = Var), ModelFit = LogMM)
    if (CalibrLines != "overall") {
      Results$EstC = predict(LogMM, newdata = PerC, re.form = ~(LP | 
                                                                  center), allow.new.levels = T, type = "response")
    }
    return(Results)
  }, simplify = F, USE.NAMES = T)
  cat("\n\nComputing overall/center-specific calibration curve(s)...\n\n")
  if (Parallel) {
    cl <- makeCluster(NrCores)
    clusterEvalQ(cl, library("REMA"))
    clusterExport(cl, c("Df", "NrImp", "NewX", "CL", "CalibrLines", 
                        "Controlglmer", "EpsGrad", "fNonConv"), envir = environment())
    MethCalc = "parSapply"
    ArgzCalc$cl = cl
    on.exit(stopCluster(cl))
  }
  else {
    MethCalc = "sapply"
  }
  ResultsCal = do.call(MethCalc, args = ArgzCalc)
  MF = lapply(ResultsCal, "[[", "ModelFit")
  Conv = unlist(lapply(MF, function(x) length(x@optinfo$conv$lme4) == 
                         0))
  if (RmNonConv) 
    ResultsCal = ResultsCal[Conv]
  ResultsB = lapply(ResultsCal, "[[", "B")
  ResultsBint = do.call("rbind", lapply(ResultsB, "[[", "Interc"))
  ResultsBslo = do.call("rbind", lapply(ResultsB, "[[", "Slope"))
  ResultsBRR = as.vector(RubinRes(list(ResultsBint, ResultsBslo))[, 
                                                                  1])
  ResultsLP = lapply(ResultsCal, "[[", "LP")
  Niter = if (CL != "PI") {
    nrow(NewX)
  }
  else {
    nrow(Df[Df$imp == 1, ])
  }
  AllPred = sapply(1:Niter, function(i) {
    do.call("rbind", lapply(ResultsLP, function(x) x[i, ]))
  }, simplify = F)
  ResultsCalRR = RubinCL(AllPred)
  if (CalibrLines != "overall") {
    ResultsC = lapply(ResultsCal, "[[", "EstC")
    ResultsC = rowMeans(do.call("cbind", ResultsC))
    ResultsC = cbind.data.frame(x = exp(PerC$LP) - 1, y = ResultsC, 
                                center = PerC$center)
    ResultsC = split(ResultsC, ResultsC$center)
  }
  if(CalibrLines != "overall"){
    par(mar = c(5,5,1,13), xpd=TRUE, pty = 's') 
  } else{
    par(pty = 's')
  }
  
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, las = las, main = title, cex.lab = 1, cex.axis = 1, ...) 
  clip(0, 1, 0, 1)
  do.call("clip", as.list(par()$usr))
  if (CalibrLines == "centers") {
    for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                              i)], lty = lty.centers[which(IncludedCenters == i)], 
                                     lwd = lwd.centers[which(IncludedCenters == i)])
    lt = c(lty.centers)
    lw.d = c(lwd.centers)
    all.col = c(col.centers)
    leg = c(as.character(IncludedCenters))
    
    lp = legendloc
    lp = list(x = lp[1], y = lp[2])
    legend("topright", leg, lty = lt, cex = 1, bty = "n", lwd = lw.d,
           col = all.col, y.intersp = y.intersp, ncol = ncol.leg, inset = c(-.58, 0), xpd=NA) 
  }
  else {
    p = if (CL == "PI") {
      ddply(Df, .(id), function(x) mean(x$p))$V1
    }
    else {
      exp(NewX$LP)-1
    }
    X = cbind(1, log(p + 1))
    FE = ResultsBRR
    y = Df$y[Df$imp == 1]
    if (!flexible) {
      OverallCal = Ilogit(X[order(p), ] %*% FE)
      p = p[order(p)]
      lines(p, OverallCal, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
      if (CL != "none") {
        lines(p, ResultsCalRR[, 2], lty = 2, col = col.overall)
        lines(p, ResultsCalRR[, 3], lty = 2, col = col.overall)
      }
    }
    else {
      Lfit = loess(y ~ p, Df)
      x = Lfit$x
      y = Lfit$fitted
      y = y[order(x)]
      x = x[order(x)]
      lines(x, y, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    lt = c(lty.overall)
    lw.d = c(lwd.overall)
    all.col = c(col.overall)
    leg = c("Overall")
    if (CalibrLines == "both") {
      for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                                i)], lty = lty.centers[which(IncludedCenters == 
                                                                                                               i)], lwd = lwd.centers[which(IncludedCenters == 
                                                                                                                                              i)])
      lt = c(lty.centers)
      lw.d = c(lwd.centers)
      all.col = c(col.centers)
      leg = c(as.character(IncludedCenters))
    }
  }
  lp = legendloc
  lp = list(x = lp[1], y = lp[2])
  if (dostats) {
    stats.2 <- paste(
      "Discrimination\n", 
      "...c-statistic: ", sprintf(paste("%.", roundstats, 
                                        "f", sep = ""), AUC$Performance[1]), " (", sprintf(paste("%.", 
                                                                                                 roundstats, "f", sep = ""), AUC$Performance[2]), 
      " to ", sprintf(paste("%.", roundstats, "f", sep = ""), 
                      AUC$Performance[3]), ")", sep = "")
    text(statloc[1], statloc[2], stats.2, pos = 4, cex = cex)
  }
  if (CalibrLines != "centers") {
    x = ddply(Df, .(id), function(x) mean(x$p))$V1
    y = Df$y[Df$imp == 1]
    bins <- seq(0, min(1, max(xlim)), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0 = table(cut(x[y == 0], bins))
    f1 = table(cut(x[y == 1], bins))
    j0 = f0 > 0
    j1 = f1 > 0
    bins0 = (bins[-101])[j0]
    bins1 = (bins[-101])[j1]
    f0 = f0[j0]
    f1 = f1[j1]
    maxf = max(f0, f1)
    f0 = (0.1 * f0)/maxf
    f1 = (0.1 * f1)/maxf
    segments(bins1, line.bins, bins1, length.seg * f1 + line.bins)
    segments(bins0, line.bins, bins0, length.seg * -f0 + 
               line.bins)
    lines(c(min(bins0, bins1) - 0.01, max(bins0, bins1) + 
              0.01), c(line.bins, line.bins))
    text(max(bins0, bins1) + dist.label, line.bins + dist.label2, 
         d1lab, cex = cex.d01)
    text(max(bins0, bins1) + dist.label, line.bins - dist.label2, 
         d0lab, cex = cex.d01)
  }
  Performance = rbind(unname(CalInterc), unname(CalSlope), 
                      unname(do.call("cbind", (AUC$Performance[, 1:3, drop = T]))))
  rownames(Performance) = c("Calibration intercept", "Calibration slope", 
                            "AUC")
  colnames(Performance) = c("Point estimate", "LCL", "UCL")
  AllResults = structure(list(Performance = Performance, AUCinfo = AUC, 
                              included = unique(Df$center), ConfLevel = 1 - alpha, 
                              ResultsBRR = ResultsBRR, ResultsCalRR = ResultsCalRR, 
                              PlotArgz = list(Mtext = list(stats = stats.2), Plot = list(x = lp, 
                                                                                         leg = leg, cex.leg = cex.leg, lwd = lw.d, col = all.col, 
                                                                                         y.intersp = y.intersp, ncol = ncol.leg)), call = Argz), 
                         class = "RE_ValProbImp")
  if (CalibrLines == "centers") {
    AllResults$PerCenter = ResultsC
  }
  else {
    AllResults$Plot = cbind.data.frame(x = p, y = OverallCal)
    if (CalibrLines == "both") 
      AllResults$PerCenter = ResultsC
  }
  return(AllResults)
  
}

## Probabilities (LR2, SRrisks, ADNEX): Pooled data
CalPool.Imp <- function(pred, outcome, imp, data,
                        title = "Calibration curve", xlab = "Predicted probability", 
                        ylab = "Observed proportion"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  LP = logit(as.numeric(eval(arguments$pred, data)))
  outcome <- eval(arguments$outcome, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, LP = LP, imp = imp, stringsAsFactors = F)
  
  NrImp = length(unique(imp))
  
  # Matrices
  FE <- matrix(ncol = 3, nrow = NrImp)
  colnames(FE) <- c("Imp", "Intercept", "LP")
  FE <- data.frame(FE)
  FE$Imp <- seq(1, NrImp, by = 1)
  
  FE.int <- matrix(ncol = 3, nrow = NrImp)
  colnames(FE.int) <- c("Imp", "Intercept", "LP")
  FE.int <- data.frame(FE.int)
  FE.int$Imp <- seq(1, NrImp, by = 1)
  
  RE <- list()
  for(i in 1:NrImp){
    RE[[i]] <- list()
  }
  

  calplotdata.avg <- as.data.frame(matrix(NA, NrImp, 2*NrImp + 1))
  colnames(calplotdata.avg) <-  c("LP", paste("obs", c(1:NrImp)), paste("var", c(1:NrImp)))
  calplotdata.avg$LP <-  seq(min(Df$LP), max(Df$LP), length.out = NrImp)
  
  
  for(i in 1:NrImp){
    cat("Imputation", i, " of", NrImp, "\n")
    cal.model <- glm(y ~ 1 + LP, data = Df[Df$imp == i, ], family = "binomial")
    FE$Intercept[i] <- coef(cal.model)[1]
    FE$LP[i] <- coef(cal.model)[2]
    FE$var_Int[i] <- vcov(cal.model)[1, 1]
    FE$var_LP[i] <- vcov(cal.model)[2, 2]
    
    cal.int <- glm(y ~ 1, offset = LP, data = Df[Df$imp == i, ], family = "binomial")
    FE.int$Intercept[i] <- coef(cal.int)[1]
    FE.int$var_Int[i] <- vcov(cal.model)
    

    # Predict probability
    calplotdata.avg[, i + 1] <- predict(cal.model, newdata = calplotdata.avg, type = "response", allow.new.levels = T)
    
    #get the variance for the confidence interval
    mm <- model.matrix( ~ calplotdata.avg$LP, calplotdata.avg)
    calplotdata.avg[, i + NrImp + 1] <-
      diag(mm %*% tcrossprod(vcov(cal.model), mm))
  }  
  
  Performance.slope <- REMA:::rubinrules(estimate = LP, variance = var_LP, data = FE, expit = T)
  Performance.int <- REMA:::rubinrules(estimate = Intercept, variance = var_Int, data = FE.int, expit = T)
  Performance <- rbind(Performance.int, Performance.slope)
  
  calplotdata.avg$p_pred <- inv.logit(calplotdata.avg$LP)
  calplotdata.avg$p_obs <- rowMeans(calplotdata.avg[ ,2:NrImp + 1], na.rm=T)
  calplotdata.avg$lp_obs <- rowMeans(calplotdata.avg[ ,2:NrImp + 1], na.rm=T)
  
  #combine between and within imputatation variance
  Bvar <- NrImp + 2
  Evar <- 2 * NrImp + 1
  Eobs <- NrImp + 1
  calplotdata.avg$var.w<-rowMeans(calplotdata.avg[ ,Bvar:Evar], na.rm=T) #gemiddelde variance binnen imputaties
  calplotdata.avg$var.b<-rowVars(as.matrix(calplotdata.avg[ ,2:Eobs]), na.rm=T) #between-imputation variance
  calplotdata.avg$var.t<-with(calplotdata.avg, var.w+(1+1/NrImp)*var.b)
  calplotdata.avg$se.t<-with(calplotdata.avg, sqrt(var.t))
  calplotdata.avg$df<-with(calplotdata.avg, (1+(NrImp/(1+NrImp))*var.w/var.b)^2*(NrImp-1))
  calplotdata.avg$t<-with(calplotdata.avg, qt(.975, df))
  calplotdata.avg$lcl<-with(calplotdata.avg, lp_obs-t*se.t)
  calplotdata.avg$ucl<-with(calplotdata.avg, lp_obs+t*se.t)
  calplotdata.avg$p_plo <- plogis(calplotdata.avg$lcl)
  calplotdata.avg$p_phi <- plogis(calplotdata.avg$ucl)
  
  
  plot(y = seq(0, 1, length.out = 100),
       x = seq(0, 1, length.out = 100),
       type = "l",
       xlab = xlab,
       ylab = ylab,
       main = title,
       col = "red",
       lwd = 1,
       xaxt = "n",
       cex.lab = 1.5,
       cex.axis = 1.5
  )
  axis(1, seq(0, 1, 0.1), labels = rep("", 11))
  lines(calplotdata.avg$p_pred,
        calplotdata.avg$p_obs,
        col = "black",
        lwd = 2)
  
  return(structure(list(Predictions = calplotdata.avg, Performance = Performance, data = Df)))
}

## Probabilities (LR2, SRrisks, ADNEX): Multiple imputed data
RE.ValProbImp <- function (p, y, center, imputation.id, patientid, data, LogCal = T, 
                           flexible = F, CL = c("none", "CI", "PI"), CalibrLines = c("overall", 
                                                                                     "centers", "both"), dostats = T, statloc = c(0, 0.85), 
                           legendloc = c(0.5, 0.27), roundstats = 2, cex = 0.75, cex.leg = 0.75, 
                           ncol.leg = 1, lty.overall = 1, lwd.overall = 2, col.overall = "red", 
                           RMprompt = F, RmNonConv = F, title = "Calibration curve", xlab = "Estimated risk of malignancy", 
                           ylab = "Observed proportion of malignancy", xlim = c(0, 1), ylim = c(0, 
                                                                                                1), d0lab = "0", d1lab = "1", cex.d01 = 0.7, dist.label = 0.04, 
                           line.bins = -0.05, dist.label2 = 0.03, las = 1, length.seg = 1, 
                           y.intersp = 1, lty.ideal = 1, col.ideal = "black", lwd.ideal = 1.75, 
                           lty.centers = NULL, lwd.centers = NULL, col.centers = NULL, 
                           Parallel = F, NrCores = detectCores() - 1, alpha = 0.05, 
                           EpsGrad = 0.001, fNonConv = c("warning", "stop"), Controlglmer = glmerControl(optimizer = "bobyqa"), 
                           AUCmeth = c("RE.auc.imp", "RE.auc.imp2"), ...) 
{
  Argz = as.list(match.call())[-1]
  p = eval(Argz$p, data)
  y = eval(Argz$y, data)
  center = eval(Argz$center, data)
  id = eval(Argz$patientid, data)
  LP = Logit(p)
  CL = match.arg(CL)
  CalibrLines = match.arg(CalibrLines)
  AUCmeth = match.arg(AUCmeth)
  if (is.factor(center)) 
    center = as.character(center)
  if (length(unique(center)) == 1) 
    stop("There is only one center, hence a random effects model should not be used here.")
  if (length(unique(center)) < 5) 
    warning("There are less than 5 centers, consider using a different method.", 
            immediate. = T)
  if (LogCal & flexible) 
    stop("LogCal and flexible cannot both be true.")
  if (missing(imputation.id)) 
    stop("Specify the variable indicating the number of the imputation!")
  imp = eval(Argz$imputation.id, data)
  if (!is.numeric(imp)) 
    imp = as.numeric(imp)
  NrImp = length(unique(imp))
  Df = data.frame(y = y, p = p, LP = Logit(p), center = center, 
                  imp = imp, id = id)
  Df = Df[with(Df, order(center, p)), ]
  RmCenter = dlply(Df[Df$imp == 1, ], .(center), function(x) if (sum(x$y == 
                                                                     1) < 10 | sum(x$y == 0) < 10) 
    unique(x$center)
    else NULL)
  RmCenter = REMA:::ListToVector(RmCenter)
  if (length(RmCenter) != 0) {
    if (RMprompt) {
      Cmmnd = readline(paste0("The center(s) ", paste0(RmCenter, 
                                                       collapse = ", "), " have less than 10 (non-)events and these will be", 
                              " removed. Do you want to continue? (y/n)   "))
      if (Cmmnd == "n") 
        stop("Function was stopped by the user.")
    }
    Df = Df[!(Df$center %in% RmCenter), ]
  }
  IncludedCenters = unique(Df$center)
  centers = IncludedCenters
  nCenters = length(IncludedCenters)
  if (!Parallel & length(unique(Df$imp)) > 50) {
    cat("\n\nParallel computing can be used to get the results faster.\n\n")
  }
  else if (Parallel & length(unique(Df$imp)) < 10) {
    cat(paste("\n\nParallel computing is only recommended when there is a large number of imputed datasets.", 
              " When used with a small number of imputed datasets, this may be slower.\n\n"))
  }
  if (CalibrLines != "overall") {
    if (all(sapply(list(lty.centers, lwd.centers, col.centers), 
                   is.null))) {
      lty.centers = rep(1, nCenters)
      lwd.centers = rep(1, nCenters)
      col.centers = seq_len(nCenters)
    }
    else if (sum(!sapply(list(lty.centers, lwd.centers, col.centers), 
                         is.null)) != 3) {
      stop("If you specify one of the arguments, please also specify the others.")
    }
    else {
      if (any(sapply(list(lty.centers, lwd.centers, col.centers), 
                     function(x) length(x) < nCenters))) 
        stop("The vector length of lty.centers, lwd.centers and col.centers is less than the number of centers.")
      FixL = function(x) x[1:nCenters]
      lty.centers = FixL(lty.centers)
      col.centers = FixL(col.centers)
      lwd.centers = FixL(lwd.centers)
    }
  }
  RubinRes <- function(x) {
    do.call("rbind", lapply(x, function(x) {
      tmp = REMA:::rubinrules(Est, Var, x, alpha = alpha)[c("Mean estimate", 
                                                            "Total variance")]
      tmp = do.call("cbind.data.frame", tmp)
      names(tmp) = c("Est", "Var")
      tmp
    }))
  }
  RubinCL <- function(x) {
    do.call("rbind", lapply(x, function(x) {
      tmp = REMA:::rubinrules(Est, Var, x, alpha = alpha, expit = T)$Expit
      names(tmp) = c("Est", "LCL", "UCL")
      tmp
    }))
  }
  NewX = data.frame(LP = seq(min(Df$LP), max(Df$LP), length = 500))
  PerC = data.frame(LP = rep(NewX$LP, length(IncludedCenters)), 
                    center = sort(rep(IncludedCenters, 500)))
  ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
    tmp = Df[Df$imp == i, ]
    LogMM = glmer(y ~ LP + (LP | center), data = tmp, family = "binomial", 
                  control = Controlglmer)
    LogMM2 = glmer(y ~ 1 + (1 | center), data = tmp, family = "binomial", 
                   offset = LP, control = Controlglmer)
    X = if (CL != "PI") {
      as.matrix(cbind(1, NewX$LP))
    } else {
      model.matrix(LogMM)[order(tmp$p), ]
    }
    FE = fixef(LogMM)
    VarFE = diag(as.matrix(vcov.merMod(LogMM)))
    Est = X %*% FE
    Var = diag(X %*% tcrossprod(as.matrix(vcov.merMod(LogMM)), 
                                X))
    FE_int = fixef(LogMM2)
    VarFE_int = diag(as.matrix(vcov.merMod(LogMM2)))
    Results = list(B = list(Interc = cbind.data.frame(Est = FE[1], 
                                                      Var = VarFE[1]), Slope = cbind.data.frame(Est = FE[2], 
                                                                                                Var = VarFE[2]), Interc2 = cbind.data.frame(Est = FE_int[1], 
                                                                                                                                            Var = VarFE_int[1])), LP = cbind.data.frame(Est = Est, 
                                                                                                                                                                                        Var = Var), ModelFit = LogMM)
    if (CalibrLines != "overall") {
      Results$EstC = predict(LogMM, newdata = PerC, re.form = ~(LP | 
                                                                  center), allow.new.levels = T, type = "response")
    }
    return(Results)
  }, simplify = F, USE.NAMES = T)
  cat("\n\nComputing overall/center-specific calibration curve(s)...\n\n")
  if (Parallel) {
    cl <- makeCluster(NrCores)
    clusterEvalQ(cl, library("REMA"))
    clusterExport(cl, c("Df", "NrImp", "NewX", "CL", "CalibrLines", 
                        "Controlglmer", "EpsGrad", "fNonConv"), envir = environment())
    MethCalc = "parSapply"
    ArgzCalc$cl = cl
    on.exit(stopCluster(cl))
  }
  else {
    MethCalc = "sapply"
  }
  ResultsCal = do.call(MethCalc, args = ArgzCalc)
  MF = lapply(ResultsCal, "[[", "ModelFit")
  Conv = unlist(lapply(MF, function(x) length(x@optinfo$conv$lme4) == 
                         0))
  if (RmNonConv) 
    ResultsCal = ResultsCal[Conv]
  ResultsB = lapply(ResultsCal, "[[", "B")
  ResultsBint = do.call("rbind", lapply(ResultsB, "[[", "Interc"))
  ResultsBint2 = do.call("rbind", lapply(ResultsB, "[[", "Interc2"))
  ResultsBslo = do.call("rbind", lapply(ResultsB, "[[", "Slope"))
  ResultsBRR = as.vector(RubinRes(list(ResultsBint, ResultsBslo))[, 
                                                                  1])
  SumSlo = RubinRes(list(ResultsBslo))
  SumSlo$LL <- SumSlo$Est - 1.96 * sqrt(SumSlo$Var)
  SumSlo$UL <- SumSlo$Est + 1.96 * sqrt(SumSlo$Var)
  
  SumInt = RubinRes(list(ResultsBint2))
  SumInt$LL <- SumInt$Est - 1.96 * sqrt(SumInt$Var)
  SumInt$UL <- SumInt$Est + 1.96 * sqrt(SumInt$Var)
  
  ResultsLP = lapply(ResultsCal, "[[", "LP")
  Niter = if (CL != "PI") {
    nrow(NewX)
  }
  else {
    nrow(Df[Df$imp == 1, ])
  }
  AllPred = sapply(1:Niter, function(i) {
    do.call("rbind", lapply(ResultsLP, function(x) x[i, ]))
  }, simplify = F)
  ResultsCalRR = RubinCL(AllPred)
  if (CalibrLines != "overall") {
    ResultsC = lapply(ResultsCal, "[[", "EstC")
    ResultsC = rowMeans(do.call("cbind", ResultsC))
    ResultsC = cbind.data.frame(x = Ilogit(PerC$LP), y = ResultsC, 
                                center = PerC$center)
    ResultsC = split(ResultsC, ResultsC$center)
  }
  
  if(CalibrLines != "overall"){
    par(mar = c(5,5,1,13), xpd=TRUE, pty = 's') 
  } else{
    par(pty = 's')
  }
  
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, las = las, main = title, cex.lab = 1, cex.axis = 1, ...) 
  clip(0, 1, 0, 1)
  abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
  do.call("clip", as.list(par()$usr))
  lt = lty.ideal
  lw.d = lwd.ideal
  all.col = col.ideal
  leg = "Ideal"
  if (CalibrLines == "centers") {
    for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                              i)], lty = lty.centers[which(IncludedCenters == i)], 
                                     lwd = lwd.centers[which(IncludedCenters == i)])
    lt = c(lt, lty.centers)
    lw.d = c(lw.d, lwd.centers)
    all.col = c(all.col, col.centers)
    leg = c(leg, as.character(IncludedCenters))
  }
  else {
    p = if (CL == "PI") {
      ddply(Df, .(id), function(x) mean(x$p))$V1
    }
    else {
      Ilogit(NewX$LP)
    }
    X = cbind(1, Logit(p))
    FE = ResultsBRR
    y = Df$y[Df$imp == 1]
    if (!flexible) {
      OverallCal = Ilogit(X[order(p), ] %*% FE)
      p = p[order(p)]
      lines(p, OverallCal, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
      if (CL != "none") {
        lines(p, ResultsCalRR[, 2], lty = 2, col = col.overall)
        lines(p, ResultsCalRR[, 3], lty = 2, col = col.overall)
      }
    }
    else {
      Lfit = loess(y ~ p, Df)
      x = Lfit$x
      y = Lfit$fitted
      y = y[order(x)]
      x = x[order(x)]
      lines(x, y, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    lt = c(lt, lty.overall)
    lw.d = c(lw.d, lwd.overall)
    all.col = c(all.col, col.overall)
    leg = c(leg, "Overall")
    if (CalibrLines == "both") {
      for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                                i)], lty = lty.centers[which(IncludedCenters == 
                                                                                                               i)], lwd = lwd.centers[which(IncludedCenters == 
                                                                                                                                              i)])
      lt = c(lt, lty.centers)
      lw.d = c(lw.d, lwd.centers)
      all.col = c(all.col, col.centers)
      leg = c(leg, as.character(IncludedCenters))
    }
  }
  lp = legendloc
  lp = list(x = lp[1], y = lp[2])
  legend("topright", leg, lty = lt, cex = 1, bty = "n", lwd = lw.d, 
         col = all.col, y.intersp = y.intersp, ncol = ncol.leg, inset = c(-.58, 0), xpd=NA) 
  if (dostats) {
    if (CalibrLines == "centers"){
      stats.2 <- paste("")
      text(statloc[1], statloc[2], stats.2, pos = 4, cex = cex)
    } else{
      stats.2 <- matrix(ncol = 2, nrow = 2)
      colnames(stats.2) <- c("", "Estimate (95% CI)")
      stats.2[1, ] <- c("Intercept", paste0(format(round(SumInt$Est, 2), nsmall = 2), " (", format(round(SumInt$LL, 2), nsmall = 2), " to ", format(round(SumInt$UL, 2), nsmall = 2), ")"))
      stats.2[2, ] <- c("Slope", paste0(format(round(SumSlo$Est, 2), nsmall = 2), " (", format(round(SumSlo$LL, 2), nsmall = 2), " to ", format(round(SumSlo$UL, 2), nsmall = 2), ")"))
      
      addtable2plot(x = statloc[1], y = statloc[2], table = stats.2, display.colnames = TRUE, cex = 0.75)
      
    }
  }
  
  if (CalibrLines != "centers") {
    x = ddply(Df, .(id), function(x) mean(x$p))$V1
    y = Df$y[Df$imp == 1]
    bins <- seq(0, min(1, max(xlim)), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0 = table(cut(x[y == 0], bins))
    f1 = table(cut(x[y == 1], bins))
    j0 = f0 > 0
    j1 = f1 > 0
    bins0 = (bins[-101])[j0]
    bins1 = (bins[-101])[j1]
    f0 = f0[j0]
    f1 = f1[j1]
    maxf = max(f0, f1)
    f0 = (0.1 * f0)/maxf
    f1 = (0.1 * f1)/maxf
    segments(bins1, line.bins, bins1, length.seg * f1 + line.bins)
    segments(bins0, line.bins, bins0, length.seg * -f0 + 
               line.bins)
    lines(c(min(bins0, bins1) - 0.01, max(bins0, bins1) + 
              0.01), c(line.bins, line.bins))
    text(max(bins0, bins1) + dist.label, line.bins + dist.label2, 
         d1lab, cex = cex.d01)
    text(max(bins0, bins1) + dist.label, line.bins - dist.label2, 
         d0lab, cex = cex.d01)
  }
  
  AllResults = structure(list( 
    included = unique(Df$center), ConfLevel = 1 - alpha, 
    ResultsBRR = ResultsBRR, ResultsCalRR = ResultsCalRR, 
    PlotArgz = list(Mtext = list(stats = stats.2), Plot = list(x = lp, 
                                                               leg = leg, cex.leg = cex.leg, lwd = lw.d, col = all.col, 
                                                               y.intersp = y.intersp, ncol = ncol.leg)), call = Argz, OneInt = ResultsBint2, OneSlo = ResultsBslo, 
    One = ResultsB, SumSlo = SumSlo, SumInt = SumInt),
    class = "RE_ValProbImp")
  if (CalibrLines == "centers") {
    AllResults$PerCenter = ResultsC
  }
  else {
    AllResults$Plot = cbind.data.frame(x = p, y = OverallCal)
    if (CalibrLines == "both") 
      AllResults$PerCenter = ResultsC
  }
  return(AllResults)
}


## Probabilities (LR2, SRrisks, ADNEX): Without multiple imputed data
RE.ValProb2 <- function (p, y, center, data, CalibrLines = c("overall", "centers", 
                                                             "both"), LogCal = T, flexible = F, dostats = T, statloc = c(0, 
                                                                                                                         0.85), legendloc = c(0.5, 0.27), MethodCL = c("profile", 
                                                                                                                                                                       "Wald", "boot"), LevelCL = 0.95, roundstats = 3, cex = 0.75, 
                         cex.leg = 0.75, ncol.leg = 1, lty.overall = 1, lwd.overall = 2, 
                         col.overall = "red", RMprompt = F, xlab = "Estimated risk", 
                         ylab = "Observed proportion", xlim = c(0, 1), ylim = c(0, 
                                                                                1), d0lab = "0", d1lab = "1", cex.d01 = 0.7, dist.label = 0.04, 
                         line.bins = -0.05, dist.label2 = 0.03, las = 1, length.seg = 1, 
                         y.intersp = 1, lty.ideal = 1, col.ideal = "black", lwd.ideal = 1.75, 
                         lty.centers = NULL, lwd.centers = NULL, col.centers = NULL, 
                         EpsGrad = 0.001, fNonConv = c("warning", "stop"), Controlglmer = glmerControl(optimizer = "bobyqa"), 
                         ...) 
{
  Argz = as.list(match.call())[-1]
  CalibrLines = match.arg(CalibrLines)
  MethCL = match.arg(MethodCL)
  p = eval(Argz$p, data)
  y = eval(Argz$y, data)
  center = eval(Argz$center, data)
  LP = Logit(p)
  if (is.factor(center)) 
    center = as.character(center)
  if (length(unique(center)) == 1) 
    stop("There is only one center, hence a random effects model should not be used here.")
  if (length(unique(center)) < 5) 
    warning("There are less than 5 centers, consider using a different method.", 
            immediate. = T)
  if (LogCal & flexible) 
    stop("LogCal and flexible cannot both be true.")
  Df = data.frame(y = y, p = p, LP = LP, center = center)
  Df = Df[with(Df, order(center, p)), ]
  RmCenter = dlply(Df, .(center), function(x) if (sum(x$y == 
                                                      1) < 10 | sum(x$y == 0) < 10) 
    unique(x$center)
    else NULL)
  RmCenter = unname(unlist(RmCenter))
  if (!is.null(RmCenter)) {
    if (RMprompt) {
      Cmmnd = readline(paste0("The center(s) ", paste0(RmCenter, 
                                                       collapse = ", "), " have less than 10 (non-)events and these will be", 
                              " removed. Do you want to continue? (y/n)   "))
      if (Cmmnd == "n") 
        stop("Function was stopped by the user.")
    }
    Df = Df[!(Df$center %in% RmCenter), ]
  }
  IncludedCenters = unique(Df$center)
  centers = IncludedCenters
  nCenters = length(IncludedCenters)
  if (CalibrLines != "overall") {
    if (all(sapply(list(lty.centers, lwd.centers, col.centers), 
                   is.null))) {
      lty.centers = rep(1, nCenters)
      lwd.centers = rep(1, nCenters)
      col.centers = seq_len(nCenters)
    }
    else if (sum(!sapply(list(lty.centers, lwd.centers, col.centers), 
                         is.null)) != 3) {
      stop("If you specify one of the arguments, please also specify the others.")
    }
    else {
      if (any(sapply(list(lty.centers, lwd.centers, col.centers), 
                     function(x) length(x) < nCenters))) 
        stop("The vector length of lty.centers, lwd.centers and col.centers is less than the number of centers.")
      FixL = function(x) x[1:nCenters]
      lty.centers = FixL(lty.centers)
      col.centers = FixL(col.centers)
      lwd.centers = FixL(lwd.centers)
    }
  }
  Df$center = factor(Df$center)
  contrasts(Df$center) = "contr.sum"
  LogMM = glmer(y ~ LP + (LP | center), data = Df, family = "binomial", 
                control = Controlglmer)
  CalSlope = fixef(LogMM)[2]
  cat("\n\nComputing confidence interval for calibration slope...\n\n")
  CalSlope = c(CalSlope, confint.merMod(LogMM, parm = "LP", 
                                        quiet = T, method = MethodCL, level = LevelCL))
  LogMM2 = glmer(y ~ 1 + (1 | center), data = Df, family = "binomial", 
                 offset = LP, control = Controlglmer)
  CalInterc = fixef(LogMM2)
  cat("\n\nComputing confidence interval for calibration intercept...\n\n")
  CalInterc = c(CalInterc, confint(LogMM2, parm = "(Intercept)", 
                                   quiet = T, method = MethodCL, level = LevelCL))
  AUC = RE.auc(p, y, center, Df, alpha = 1 - LevelCL)
  
  if(CalibrLines != "overall"){
    par(mar = c(5,5,1,13), xpd=TRUE, pty = 's') 
  } else{
    par(pty = 's')
  }
  
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, las = las, ...)
  clip(0, 1, 0, 1)
  abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
  do.call("clip", as.list(par()$usr))
  lt = lty.ideal
  lw.d = lwd.ideal
  all.col = col.ideal
  leg = "Ideal"
  NewX = data.frame(LP = seq(min(Df$LP), max(Df$LP), length = 500))
  p = Ilogit(NewX$LP)
  X = cbind(1, Logit(p))
  FE = fixef(LogMM)
  if (CalibrLines != "overall") {
    PerC = data.frame(LP = rep(NewX$LP, length(IncludedCenters)), 
                      center = sort(rep(IncludedCenters, 500)))
    EstC = predict(LogMM, newdata = PerC, re.form = ~(LP | 
                                                        center), allow.new.levels = T, type = "response")
    ResultsC = cbind.data.frame(x = Ilogit(PerC$LP), y = EstC, 
                                center = PerC$center)
    ResultsC = split(ResultsC, ResultsC$center)
  }
  if (CalibrLines == "centers") {
    for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                              i)], lty = lty.centers[which(IncludedCenters == i)], 
                                     lwd = lwd.centers[which(IncludedCenters == i)])
    lt = c(lt, lty.centers)
    lw.d = c(lw.d, lwd.centers)
    all.col = c(all.col, col.centers)
    leg = c(leg, as.character(IncludedCenters))
  }
  else {
    p = Ilogit(NewX$LP)
    X = cbind(1, Logit(p))
    FE = fixef(LogMM)
    if (!flexible) {
      OverallCal = Ilogit(X[order(p), ] %*% FE)
      p = p[order(p)]
      lines(p, OverallCal, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    else {
      Lfit = loess(y ~ p, Df)
      x = Lfit$x
      y = Lfit$fitted
      y = y[order(x)]
      x = x[order(x)]
      lines(x, y, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    lt = c(lt, lty.overall)
    lw.d = c(lw.d, lwd.overall)
    all.col = c(all.col, col.overall)
    leg = c(leg, "Overall")
    if (CalibrLines == "both") {
      for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                                i)], lty = lty.centers[which(IncludedCenters == 
                                                                                                               i)], lwd = lwd.centers[which(IncludedCenters == 
                                                                                                                                              i)])
      lt = c(lt, lty.centers)
      lw.d = c(lw.d, lwd.centers)
      all.col = c(all.col, col.centers)
      leg = c(leg, as.character(IncludedCenters))
    }
  }
  lp = legendloc
  lp = list(x = lp[1], y = lp[2])
  legend("topright", leg, lty = lt, cex = cex.leg, bty = "n", lwd = lw.d, 
         col = all.col, y.intersp = y.intersp, ncol = ncol.leg, inset = c(-.68, 0), xpd=NA) # Oorspronkelijk: c(-.39, 0)
  if (dostats) {
    if (CalibrLines == "centers"){
      stats.2 <- paste("")
      text(statloc[1], statloc[2], stats.2, pos = 4, cex = 0.75)
    } else{
      stats.2 <- matrix(ncol = 2, nrow = 2)
      colnames(stats.2) <- c("", "Estimate (95% CI)")
      stats.2[1, ] <- c("Intercept", paste0(format(round(CalInterc[1], 2), nsmall = 2), " (", format(round(CalInterc[2], 2), nsmall = 2), " to ", format(round(CalInterc[3], 2), nsmall = 2), ")"))
      stats.2[2, ] <- c("Slope", paste0(format(round(CalSlope[1], 2), nsmall = 2), " (", format(round(CalSlope[2], 2), nsmall = 2), " to ", format(round(CalSlope[3], 2), nsmall = 2), ")"))
      
      addtable2plot(x = statloc[1], y = statloc[2], table = stats.2, display.colnames = TRUE, cex = 0.75)
      
    }
  }
  
  Performance = rbind(CalInterc, CalSlope, AUC$Performance)
  rownames(Performance) = c("Calibration intercept", "Calibration slope", 
                            "AUC")
  colnames(Performance) = c("Point estimate", "LCL", "UCL")
  return(structure(list(Performance = Performance, AUCinfo = AUC, 
                        included = unique(Df$center), FitCalIntercept = LogMM2, 
                        FitCalSlope = LogMM, ConfLevel = LevelCL), class = "RE_ValProb"))
}


########################################
### 4. Function for data for WinBugs ###
########################################

## Without multiple imputation
DataWinBugs <- function(pred, outcome, center, data, 
                        sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  
  for(i in 1:length(sequence)){
    threshold <- sequence[i]
    
    Confusion <- matrix(nrow = length(centers), ncol = 8)
    Confusion <- data.frame(Confusion)
    colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
    Confusion$CutOff <- threshold
    
    for(j in seq_along(centers)){
      
      Confusion$Center[j] <- centers[j]
      
      CM <- confusion.matrix(obs = Df$y[Df$center == centers[j]], pred = Df$p[Df$center == centers[j]], threshold = threshold)
      Confusion$TN[j] <- CM[1,1]
      Confusion$TP[j] <- CM[2,2]
      Confusion$FP[j] <- CM[2,1]
      Confusion$FN[j] <- CM[1,2]
      
      Confusion$cases <- Confusion$TP + Confusion$FN
      Confusion$controls <- Confusion$TN + Confusion$FP
      
      Confusion$n <- Confusion$cases + Confusion$controls
      Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
      
    }
    
    ConfusionList[[i]] <- Confusion
    
    
  }
  return(structure(list(Results = ConfusionList)))
}


## With multiple imputation
DataWinBugs.imp <- function(pred, outcome, center, data, imp,
                            sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  NRimp <- length(unique(Df$imp))
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  ConfusionImp <- list()
  for(i in 1:NRimp){
    ConfusionImp[[i]] <- list()
  }
  
  for(k in 1:NRimp){
    cat("Imputation", k, "of", NRimp)
    for(i in 1:length(sequence)){
      threshold <- sequence[i]
      
      Confusion <- matrix(nrow = length(centers), ncol = 8)
      Confusion <- data.frame(Confusion)
      colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
      Confusion$CutOff <- threshold
      
      for(j in seq_along(centers)){
        
        Confusion$Center[j] <- centers[j]
        
        CM <- confusion.matrix(obs = Df$y[Df$center == centers[j] & Df$imp == k], pred = Df$p[Df$center == centers[j] & Df$imp == k], threshold = threshold)
        Confusion$TN[j] <- CM[1,1]
        Confusion$TP[j] <- CM[2,2]
        Confusion$FP[j] <- CM[2,1]
        Confusion$FN[j] <- CM[1,2]
        
        Confusion$cases <- Confusion$TP + Confusion$FN
        Confusion$controls <- Confusion$TN + Confusion$FP
        
        Confusion$n <- Confusion$cases + Confusion$controls
        Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
        
      }
      
      ConfusionList[[i]] <- Confusion
    }
    cat("rbindlist toepassen")
    ConfusionImp[[k]] <- rbindlist(ConfusionList, fill = TRUE)
  }
  
  ConfusionLong <- rbindlist(ConfusionImp, fill = TRUE)
  ConfusionSum <- summaryBy(cbind(TP, TN, cases, controls, FP, FN, NB) ~ cbind(CutOff, Center), data = ConfusionLong, FUN = mean)
  
  return(structure(list(Results = ConfusionSum)))
}

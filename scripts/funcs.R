###############
## Functions ##
###############

## Source parameters for choice of dataset
source("scripts/pars.R")

## Function to calculate age midpoint given a lower & upper age bound
get_midpoint <- function (age_min, age_max) { 
  
  midpoint <- age_min + (age_max - age_min) / 2
  
  return(midpoint)
  
}


## Function to calculate the joint binomial likelihood 
if(P$data_choice == "Vieira") {
  
  # For fitting to Vieira et al. (2015) data
  loglik <- function (kIgG, nIgG, prevIgG, kTgERP, nTgERP, prevTgERP) {
    
    dbinom(kIgG, nIgG, prevIgG, log=T) + dbinom(kTgERP, nTgERP, prevTgERP, log=T)
    
  }
  
} else if (P$data_choice == "Mangiavacchi") {
  
  # For fitting to Mangiavacchi et al. (2016) data
  loglik <- function (kIgG, nIgG, prevIgG, kTgERPsd, nTgERPsd, prevTgERPsd, kTgERProc, nTgERProc, prevTgERProc) {
    
    dbinom(kIgG, nIgG, prevIgG, log=T) + dbinom(kTgERPsd, nTgERPsd, prevTgERPsd, log=T) + dbinom(kTgERProc, nTgERProc, prevTgERProc, log=T)
    
  }
  
}


## Likelihood wrapper function
loglik.wrap <- function (par, data) {
  
  P$log.lambda0  <- par[1]
  P$log.lambda1  <- par[2]
  P$log.gradient <- par[3]
  P$log.tgerpR   <- par[4]
  
  # For fitting to Vieira et al. (2015) data
  if (P$data_choice == "Vieira") {
    
    # For Vieira data
    sol       <- ode(y = state, times = age, parms = P,  func = ab_mod)
    mod_igg   <- sol[2:6,"op_iggP"]
    mod_tgerp <- sol[2:6,"op_tgerpP_sd"]
    
    logliks <- loglik(kIgG = df$k_igg, kTgERP= df$k_tgerp_sd, nIgG = df$n, nTgERP = df$n,
                      prevIgG =  mod_igg, prevTgERP = mod_tgerp)
    
    # For fitting to Mangiavacchi et al. (2016) data
  } else if (P$data_choice == "Mangiavacchi") {
    
    sol           <- ode(y = state, times = age, parms = P,  func = ab_mod)
    mod_igg       <- sol[2:5,"op_iggP"]
    mod_tgerp_roc <- sol[2:5,"op_tgerpP_roc"]
    mod_tgerp_sd  <- sol[2:5,"op_tgerpP_sd"]
    
    logliks <- loglik(kIgG = df$k_igg, kTgERProc = df$k_tgerp_roc, kTgERPsd = df$k_tgerp_sd,
                      nIgG = df$n, nTgERProc = df$n, nTgERPsd = df$n,
                      prevIgG = mod_igg, prevTgERProc = mod_tgerp_roc, prevTgERPsd = mod_tgerp_sd)
  }

return(sum(-logliks))
  
}


## Maximum likelihood fitting function
fitmod <- function () {
  
  # Set parameter names for ML fitting
  parnames(loglik.wrap) <- c("log.lambda0","log.lambda1", "log.gradient", "log.tgerpR")
  
  # Model 1: 4 unfixed parameters: lambda0, lambda1, gradient and tgerpR (fixed = none)
  mod1 <- mle2(minuslogl=loglik.wrap, 
               start=c(log.lambda0=log(0.02), log.lambda1=log(0.02), 
                       log.gradient=log(0.05), log.tgerpR=log(0.05)), 
               method="L-BFGS-B", lower=c(-Inf, -Inf, -Inf, -Inf), 
               upper=c(Inf, Inf, Inf, Inf), data=list(df))
  
  # Model 2: 3 unfixed parameters: lambda0, lambda1 and tgerpR (fixed = gradient)
  mod2 <- mle2(minuslogl=loglik.wrap, 
               start=c(log.lambda0=log(0.02), log.lambda1=log(0.02), 
                       log.tgerpR=log(0.05)), 
               fixed=list(log.gradient=-Inf), method="L-BFGS-B", 
               lower=c(-Inf, -Inf, -Inf), upper=c(Inf, Inf, Inf), data=list(df))
  
  # Model 3: 2 unfixed parameters: lambda1 and tgerpR (fixed = lambda0 and gradient)
  mod3 <- mle2(minuslogl=loglik.wrap, 
               start=c(log.lambda1=log(0.02), log.tgerpR=log(0.05)), 
               fixed=list(log.lambda0=-Inf, log.gradient=-Inf), method="L-BFGS-B",  
               lower=c(-Inf, -Inf), upper=c(Inf, Inf), data=list(df))
  
  # Model 4: 2 unfixed parameters: lambda0 and tgerpR (fixed = lambda1 and gradient)
  mod4 <- mle2(minuslogl=loglik.wrap, start=c(log.lambda0=log(0.02), log.tgerpR=log(0.05)), 
               fixed=list(log.lambda1=-Inf, log.gradient=-Inf), method="L-BFGS-B", 
               lower=c(-Inf, -Inf), upper=c(Inf, Inf), data=list(df))
  
  named_list <- list(mod1=mod1, mod2=mod2, mod3=mod3, mod4=mod4)
  return(named_list)
  
}


## Likelihood ratio test for nested models (where diff in no. parameters = 1)
lrtest <- function (loglik_nested, loglik_complex) { 
  
  teststat <- -2 * (as.numeric(loglik_nested) - as.numeric(loglik_complex))
  pvalue <- pchisq(teststat, df=1, lower.tail=F)
  
  return(pvalue)
  
}


## Function to select best-fitting model based on P values (numeric vector of length 3) derived from LRT
select_bestfit <- function (pvals) {
  
  # If only one P value > 0.05, select model following (e.g. if pvals[1] is non-significant, select model 2)
  if (length(which(pvals > 0.05)) == 1) {
    
    if (which(pvals > 0.05) == 1){
      bestfit <- fit$mod2
    } else if (which(pvals > 0.05) == 2) { 
      bestfit <- fit$mod3
    } else if (which(pvals > 0.05) == 3) { 
      bestfit <- fit$mod4
    }
    
    # But if >1 P value > 0.05, select model following FINAL non-significant P value
  } else if (length(which(pvals > 0.05)) > 1) {
    
    if (pvals[1] > 0.05 & pvals[2] > 0.05) { 
      bestfit <- fit$mod3
    } else if (pvals[3] > 0.05 ) { 
      bestfit <- fit$mod4
    }
    
  }
  return(bestfit)
  
}


## Function to extract parameter point estimates and 95% CIs from best-fiting model object
getfit <- function (mod) {
  
  # Define vectors
  lambda0_lower <- lambda0_upper <- lambda1_lower <- lambda1_upper <- 
    gradient_lower <- gradient_upper <- tgerp_lower <- tgerp_upper <- vector("numeric", 1)
  
  # Extract point estimates of parameters
  lambda0  <- exp(coef(mod)["log.lambda0"])
  lambda1  <- exp(coef(mod)["log.lambda1"])
  gradient <- exp(coef(mod)["log.gradient"])
  tgerp    <- 1/exp(coef(mod)["log.tgerpR"])
  
  # Extract log likelihood
  loglik <- logLik(mod)[1]
  
  # Extract 95% CIs of parameter estimates
  ci1 <- profile(mod)
  
  if (lambda0 != 0) {
    lambda0_lower  <- exp(confint(ci1))["log.lambda0",1]
    lambda0_upper  <- exp(confint(ci1))["log.lambda0",2]
  }
  
  if (lambda1 != 0) {
    lambda1_lower  <- exp(confint(ci1))["log.lambda1",1]
    lambda1_upper  <- exp(confint(ci1))["log.lambda1",2]
  }
  
  if (gradient != 0) {
    gradient_lower <- exp(confint(ci1))["log.gradient",1]
    gradient_upper <- exp(confint(ci1))["log.gradient",2]
  }
  
  tgerp_lower      <- 1/exp(confint(ci1))["log.tgerpR",2]
  tgerp_upper      <- 1/exp(confint(ci1))["log.tgerpR",1]
  
  estimates <- list(lambda0=lambda0, lambda0_lower=lambda0_lower, lambda0_upper=lambda0_upper, 
                    lambda1=lambda1, lambda1_lower=lambda1_lower, lambda1_upper=lambda1_upper, 
                    gradient=gradient, gradient_lower=gradient_lower, gradient_upper=gradient_upper, 
                    tgerp=tgerp, tgerp_lower=tgerp_lower, tgerp_upper=tgerp_upper, loglik=loglik)
  
  return(estimates)
  
}


## Data cleaning function used after sensitivity analysis
clean_data <- function (data) {
  
  data <- data %>%
    na_if(., 0) %>%                 #convert 0 to NA
    select_if(~!all(is.na(.))) %>%  #remove NA-only columns
    round(., 3) 
  
  return(data)
  
}

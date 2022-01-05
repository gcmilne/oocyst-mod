#######################
## Fit model to data ##
#######################

## Load packages
library(deSolve)
library(bbmle)
library(MASS)

## Load scripts
source("scripts/pars.R")
source("scripts/data.R")
source("scripts/model.R")

## Maximum likelihood estimation of parameter values
fit <- fitmod()

## Compare model fits
logliks <- as.numeric(lapply(fit, logLik)) #store log likelihood

pvals <- vector("numeric", length=3)  #vector to store P values from lrtest

# mod2 vs. mod1
pvals[1] <- lrtest(loglik_nested=logliks[2], loglik_complex=logliks[1])

# mod3 vs. mod2
pvals[2] <- lrtest(loglik_nested=logliks[3], loglik_complex=logliks[2])

# mod4 vs. mod3
pvals[3] <- lrtest(loglik_nested=logliks[4], loglik_complex=logliks[3])

## Calculate best-fit model
bestfit <- select_bestfit(pvals)

## Estimate uncertainty in best-fitting model estimates
mod_coef <- coef(bestfit)  #extract parameter coefficients
mod_vcov <- vcov(bestfit)  #calculate variance-covariance matrix

# Set fixed parameter values
P[names(mod_coef)] <- mod_coef
P$rho              <- 1
P$se_tgerp_sd      <- 0.80
P$sp_tgerp_sd      <- 0.95

# Create sets of parameter values
set.seed(1001)
nsim     <- 1000
mod_coef <- coef(bestfit)[which(coef(bestfit) != -Inf)]  #remove pars that equal 0
par_set  <- mvrnorm(nsim, mu=mod_coef, Sigma=mod_vcov)
store    <- vector("list", nsim)

# Sample parameter values from sets, run model & store output
if (P$data_choice == "Vieira") { 
  
  for (i in 1:nrow(par_set)) { 
    P$log.lambda0 <- par_set[i,1]
    P$log.tgerpR  <- par_set[i,2]
    sol <- ode(y=state, times=age, parms=P, func=ab_mod)
    age            <- sol[,"time"         ]
    mod_igg        <- sol[,"op_iggP"      ]
    mod_tgerp_sd   <- sol[,"op_tgerpP_sd" ]
    store[[i]] <- data.frame(age=age, mod_igg, mod_tgerp_sd, sim=i) 
  }
  
} else if (P$data_choice == "Mangiavacchi") {
  
  for (i in 1:nrow(par_set)) { 
    P$log.lambda0 <- par_set[i,1]
    P$log.lambda1 <- par_set[i,2]
    P$log.tgerpR  <- par_set[i,3]
    sol <- ode(y=state, times=age, parms=P, func=ab_mod)
    age            <- sol[,"time"         ]
    mod_igg        <- sol[,"op_iggP"      ]
    mod_tgerp_roc  <- sol[,"op_tgerpP_roc"]
    mod_tgerp_sd   <- sol[,"op_tgerpP_sd" ]
    store[[i]] <- data.frame(age=age, mod_igg, mod_tgerp_roc, mod_tgerp_sd, sim=i) 
  }
  
}

out <- do.call(rbind, store)
out <- out[out$age != 0,]  #remove rows where age=0

# Extract means, lower and upper 95% CIs
res_ig   <- aggregate(out$mod_igg,      list(out$age), FUN=quantile, probs=c(0.025, 0.50, 0.975)) #IgG
res_tgsd <- aggregate(out$mod_tgerp_sd, list(out$age), FUN=quantile, probs=c(0.025, 0.50, 0.975)) #TgERP SD

if (P$data_choice == "Mangiavacchi") {
  res_tgroc <- aggregate(out$mod_tgerp_roc, list(out$age), FUN=quantile, probs=c(0.025, 0.50, 0.975)) #TgERP ROC
}

# Create dataframe for plotting
if (P$data_choice == "Vieira") {
  v_fit <- data.frame(
    age        = df$age,
    ig_mean    = res_ig$x[,2], 
    ig_lower   = res_ig$x[,1], 
    ig_upper   = res_ig$x[,3],
    tgsd_mean  = res_tgsd$x[,2],
    tgsd_lower = res_tgsd$x[,1],
    tgsd_upper = res_tgsd$x[,3]
  )
  rm(res_ig, res_tgsd) #remove objects from environment
  
} else if (P$data_choice == "Mangiavacchi") {
  m_fit <- data.frame(
    age         = df$age,
    ig_mean     = res_ig$x[,2], 
    ig_lower    = res_ig$x[,1], 
    ig_upper    = res_ig$x[,3],
    tgsd_mean   = res_tgsd$x[,2],
    tgsd_lower  = res_tgsd$x[,1],
    tgsd_upper  = res_tgsd$x[,3], 
    tgroc_mean  = res_tgroc$x[,2],
    tgroc_lower = res_tgroc$x[,1],
    tgroc_upper = res_tgroc$x[,3]
  )
  rm(res_ig, res_tgsd, res_tgroc) #remove objects from environment
}

# Save dataframe
if (P$data_choice == "Vieira") {
  saveRDS(v_fit, paste("mod_output/", P$data_choice, "_modfit.RDS", sep=""))
  
} else if (P$data_choice == "Mangiavacchi") {
  saveRDS(m_fit, paste("mod_output/", P$data_choice, "_modfit.RDS", sep=""))
}

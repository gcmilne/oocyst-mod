##########################
## Sensitivity analyses ##
##########################

## Load packages
library(deSolve)
library(bbmle)
library(MASS)
library(purrr)
library(dplyr)

## Load scripts
source("scripts/pars.R")
source("scripts/data.R")
source("scripts/model.R")

## Sensitivity analysis 1: Fixed proportion of oocyst infections, varying diagnostic se & sp ##

# Create set up of sensitivity & specificity values
par_sets    <- setNames(data.frame(matrix(ncol=2, nrow=9)), c("se", "sp"))
par_sets$se <- rep(seq(0.80, 0.90, 0.05), each=3)
par_sets$sp <- rep(seq(0.85, 0.95, 0.05), times=3)

# Set fixed parameter values
P$rho <- 1

# To store parameter estimates & 95% CIs
sens1 <- setNames(data.frame(matrix(ncol=13, nrow=nrow(par_sets))), 
                  c("lambda0", "lambda0_lower", "lambda0_upper", 
                    "lambda1", "lambda1_lower", "lambda1_upper", 
                    "gradient", "gradient_lower", "gradient_upper", 
                    "tgerp", "tgerp_lower", "tgerp_upper", "loglik"))

# To store all model log likelihoods
logliks1 <- vector("list", length=nrow(par_sets))

# To store P values for likelihood ratio test model comparisons
pvals <- vector("numeric", length=3)

# Do sensitivity analysis
for(i in 1:nrow(par_sets)) {
  
  # Change par set
  P$se_tgerp_sd <- par_sets$se[i]
  P$sp_tgerp_sd <- par_sets$sp[i]
  
  ## Perform fitting
  fit <- fitmod()
  
  ## Compare model fits
  # Store log likelihood
  logliks1[[i]] <- as.numeric(lapply(fit, logLik))
  
  # mod2 vs. mod1
  pvals[1] <- lrtest(loglik_nested=logliks1[[i]][2], loglik_complex=logliks1[[i]][1])
  
  # mod3 vs. mod2
  pvals[2] <- lrtest(loglik_nested=logliks1[[i]][3], loglik_complex=logliks1[[i]][2])
  
  # mod4 vs. mod3
  pvals[3] <- lrtest(loglik_nested=logliks1[[i]][4], loglik_complex=logliks1[[i]][3])
  
  ## Calculate best-fit model
  bestfit <- select_bestfit(pvals)
  
  ## Extract parameter point estimates & 95% CIs, & model log likelihood
  est <- getfit(bestfit)
  
  # lambda0
  sens1$lambda0[i]        <- est$lambda0
  sens1$lambda0_lower[i]  <- est$lambda0_lower
  sens1$lambda0_upper[i]  <- est$lambda0_upper
  
  # lambda1
  sens1$lambda1[i]        <- est$lambda1
  sens1$lambda1_lower[i]  <- est$lambda1_lower
  sens1$lambda1_upper[i]  <- est$lambda1_upper
  
  # gradient
  sens1$gradient[i]       <- est$gradient
  sens1$gradient_lower[i] <- est$gradient_lower
  sens1$gradient_upper[i] <- est$gradient_upper
  
  # TgERP reversion, given as a duration (1/rate)
  sens1$tgerp[i]          <- est$tgerp
  sens1$tgerp_lower[i]    <- est$tgerp_lower
  sens1$tgerp_upper[i]    <- est$tgerp_upper
  
  # log likelihood
  sens1$loglik[i]         <- est$loglik  #log likelihood
  
}


## Clean output data
sens1 <- clean_data(sens1)

logliks1 <- logliks1 %>% 
  map(., round, 2)  #round values

# Save the results
saveRDS(sens1, paste("mod_output/", P$data_choice, "_sens1.RDS", sep=""))
saveRDS(logliks1, paste("mod_output/", P$data_choice, "_sens1_logliks.RDS", sep=""))


## Sensitivity analysis 2: Fixed diagnostic se & sp, varying proportion of oocyst infections ##
# Set fixed parameter values
P$se_tgerp_sd <- 0.80
P$sp_tgerp_sd <- 0.95

# Create set up of sensitivity & specificity values
rho_vec <- seq(0.5, 1, 0.1)

# To store parameter estimates & 95% CIs
sens2 <- setNames(data.frame(matrix(ncol=13, nrow=length(rho_vec))), 
                  c("lambda0", "lambda0_lower", "lambda0_upper", 
                    "lambda1", "lambda1_lower", "lambda1_upper", 
                    "gradient", "gradient_lower", "gradient_upper", 
                    "tgerp", "tgerp_lower", "tgerp_upper", "loglik"))

# To store all model log likelihoods
logliks2 <- vector("list", length=length(rho_vec))

# To store P values for likelihood ratio test model comparisons
pvals <- vector("numeric", length=3)

## Do sensitivity analysis
for (i in 1:length(rho_vec)) {
  
  ## Change par set
  P$rho <- rho_vec[i]
  
  ## Perform fitting
  fit <- fitmod()
  
  ## Compare model fits
  # Store log likelihood
  logliks2[[i]] <- as.numeric(lapply(fit, logLik))
  
  # mod2 vs. mod1
  pvals[1] <- lrtest(loglik_nested=logliks2[[i]][2], loglik_complex=logliks2[[i]][1])
  
  # mod3 vs. mod2
  pvals[2] <- lrtest(loglik_nested=logliks2[[i]][3], loglik_complex=logliks2[[i]][2])
  
  # mod4 vs. mod3
  pvals[3] <- lrtest(loglik_nested=logliks2[[i]][4], loglik_complex=logliks2[[i]][3])
  
  ## Calculate best-fit model
  bestfit <- select_bestfit(pvals)
  
  ## Extract parameter point estimates & 95% CIs, & model log likelihood
  est <- getfit(bestfit)
  
  # lambda0
  sens2$lambda0[i]        <- est$lambda0
  sens2$lambda0_lower[i]  <- est$lambda0_lower
  sens2$lambda0_upper[i]  <- est$lambda0_upper
  
  # lambda1
  sens2$lambda1[i]        <- est$lambda1
  sens2$lambda1_lower[i]  <- est$lambda1_lower
  sens2$lambda1_upper[i]  <- est$lambda1_upper
  
  # gradient
  sens2$gradient[i]       <- est$gradient
  sens2$gradient_lower[i] <- est$gradient_lower
  sens2$gradient_upper[i] <- est$gradient_upper
  
  # TgERP reversion, given as a duration (1/rate)
  sens2$tgerp[i]          <- est$tgerp
  sens2$tgerp_lower[i]    <- est$tgerp_lower
  sens2$tgerp_upper[i]    <- est$tgerp_upper
  
  # log likelihood
  sens2$loglik[i]         <- est$loglik  #log likelihood
  
}

## Clean output data
sens2 <- clean_data(sens2)

logliks2 <- logliks2 %>% 
  map(., round, 2)  #round values

# Save the results
saveRDS(sens2, paste("mod_output/", P$data_choice, "_sens2.RDS", sep=""))
saveRDS(logliks2, paste("mod_output/", P$data_choice, "_sens2_logliks.RDS", sep=""))

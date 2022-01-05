##########
## Data ##
##########

## Load packages
library("binom")

## Load scripts
source("scripts/funcs.R")

## Create seroprevalence datasets
if (P$data_choice == "Vieira") {  #Vieira et al. (2015)
  
  df <- setNames(data.frame(matrix(nrow=5, ncol=4)), c("age","k_igg", "k_tgerp_sd", "n"))
  df$age        <- c(get_midpoint(8, 19), get_midpoint(20, 29), get_midpoint(30, 39), #age midpoints
                     get_midpoint(40, 49), get_midpoint(50, 59)) 
  df$k_igg      <- c(11, 18, 19, 11, 22) #no. non-stage-specific IgG-positives 
  df$k_tgerp_sd <- c(9, 15, 11, 7, 12)   #no. anti-TgERP IgG-positives
  df$n          <- c(23, 27, 21, 15, 22) #sample size (n)

} else if (P$data_choice == "Mangiavacchi") { #Mangiavacchi et al. (2016)
  
  df <- setNames(data.frame(matrix(nrow=4, ncol=5)),  c("age","k_igg", "k_tgerp_sd", "k_tgerp_roc", "n"))
  df$age         <- c(get_midpoint(0, 7), get_midpoint(8, 14),  #age midpoints
                      get_midpoint(15, 21), get_midpoint(22, 28)) 
  df$k_igg       <- c(13, 59, 33, 45)   #no. non-stage-specific IgG-positives 
  df$k_tgerp_sd  <- c(9, 20, 13, 31)    #no. anti-TgERP IgG-positives (SD)
  df$k_tgerp_roc <- c(36, 75, 34, 45)   #no. anti-TgERP IgG-positives (ROC)
  df$n           <- c(116, 156, 48, 60) #sample size (n)
  
}

## Calculate 95% confidence intervals for seroprevalence data
if (P$data_choice == "Vieira") {
  cis <- binom.confint(x=c(df$k_igg, df$k_tgerp_sd), n=rep(df$n, 2), conf.level=0.95, methods="ac")

} else if (P$data_choice == "Mangiavacchi") { 
  cis <- binom.confint(x=c(df$k_igg, df$k_tgerp_sd, df$k_tgerp_roc), n=rep(df$n, 3), conf.level=0.95, methods="ac")
}

cis$upper[which(cis$upper > 1)] <- 1 #limit upper bound to 1

df$ig_mean    <- df$k_igg/df$n
df$ig_lower   <- cis$lower[1:length(df$k_igg)]
df$ig_upper   <- cis$upper[1:length(df$k_igg)]
df$tgsd_mean  <- df$k_tgerp_sd/df$n
df$tgsd_lower <- cis$lower[(length(df$k_igg)+1):(length(df$k_igg)*2)]
df$tgsd_upper <- cis$upper[(length(df$k_igg)+1):(length(df$k_igg)*2)]

if (P$data_choice == "Mangiavacchi") { 
  df$tgroc_mean  <- df$k_tgerp_roc/df$n
  df$tgroc_lower <- cis$lower[(length(df$k_igg)*2+1):(length(df$k_igg)*3)]
  df$tgroc_upper <- cis$upper[(length(df$k_igg)*2+1):(length(df$k_igg)*3)]
}

## Define age midpoints for model
age <- c(0, df$age)

################
## Parameters ##
################

## Toggle to choose dataset to fit
data_choice_options <- c("Vieira", "Mangiavacchi")

## Define parameter vector
P     <- list(
  
  data_choice = data_choice_options[1], #change from [1] to [2] to change dataset
  
  ## Fixed-rate parameters: ##
  iggR = 0,   #reversion rate of non-stage-specific IgG
  se_igg       = mean(c(0.989, 0.956, 1.00, 0.994)),      #sensitivity of conventional ELISA
  sp_igg       = mean(c(0.983, 0.987, 0.985, 0.994)),     #specificity of conventional ELISA
  se_tgerp_roc = mean(c(0.9231, 0.9153, 0.9091, 0.9111)), #sensitivity of TgERP (ROC) ELISA
  sp_tgerp_roc = mean(c(0.7670, 0.7835, 0.7333, 0.7333)), #specificity of TgERP (ROC) ELISA
  
  ## Parameters to perform sensitivity analyses on: ##
  se_tgerp_sd = 0.80,  #sensitivity of TgERP (SD) ELISA - to be varied in sensitivity analysis
  sp_tgerp_sd = 0.95,  #specificity of TgERP (SD) ELISA - to be varied in sensitivity analysis
  rho         = 1.0,   #proportion of infections that are oocyst-acquired
  
  ## Parameters to be fit to data: ##
  log.tgerpR   = log(0.129157998),  #reversion rate of anti-TgERP IgG
  log.lambda0  = log(0.01),         #baseline FoI 
  log.lambda1  = log(0.001),        #interaction of FoI with age
  log.gradient = log(0.06)          #determines shape of force of infection (FoI) over age
  
)

# Define state variables & initial conditions
state <- c(tp_iggP   = 0,  #true prevalence (tp) of non-stage-specific IgG-positives
           tp_iggN   = 1,  #true prevalence (tp) of non-stage-specificIgG-negatives
           tp_tgerpP = 0,  #true prevalence (tp) of anti-TgERP IgG-positives
           tp_tgerpN = 1)  #true prevalence (tp) of anti-TgERP IgG-negatives

rm(data_choice_options)  #clear from environment

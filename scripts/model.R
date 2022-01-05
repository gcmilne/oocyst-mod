#########################
## Serocatalytic model ##
#########################

# Tracks change in seroprevalence of non-stage-specific IgG & anti-TgERP IgG over age

ab_mod = function(age, state, P)
{
  with(as.list(c(age, state, P)),
       {
         
         # Take exponents of log-transformed parameters
         tgerpR   <- exp(P$log.tgerpR)
         gradient <- exp(P$log.gradient)
         lambda0  <- exp(P$log.lambda0)
         lambda1  <- exp(P$log.lambda1)
         
         # Force of infection (FoI) functional form
         foi <- lambda0 + lambda1 * (age * exp(-gradient*age))
         
         # Ordinary differential equations (ODEs)
         dtp_iggP   <- (       foi  *   tp_iggN ) -   iggR *   tp_iggP
         dtp_iggN   <- (      -foi  *   tp_iggN ) +   iggR *   tp_iggP
         dtp_tgerpP <- (rho *  foi) * tp_tgerpN   - tgerpR * tp_tgerpP
         dtp_tgerpN <- (rho * -foi) * tp_tgerpN   + tgerpR * tp_tgerpP
         
         # Adjust true seroprevalence to observed seroprevalence (Diggle 2011)
         op_iggP       <-  tp_iggP   * (se_igg       + sp_igg       -1) + (1 - sp_igg       )
         op_tgerpP_roc <-  tp_tgerpP * (se_tgerp_roc + sp_tgerp_roc -1) + (1 - sp_tgerp_roc )
         op_tgerpP_sd  <-  tp_tgerpP * (se_tgerp_sd  + sp_tgerp_sd  -1) + (1 - sp_tgerp_sd  )
         
         # Return output
         return(list(c(dtp_iggP, dtp_iggN, dtp_tgerpP, dtp_tgerpN), 
                     op_iggP=op_iggP, op_tgerpP_roc=op_tgerpP_roc, op_tgerpP_sd=op_tgerpP_sd, foi=foi))
       })
}

##########################
## Plot model estimates ##
##########################

## Load packages
library(ggplot2)
library(ggpubr)

## Load scripts
source("scripts/pars.R")
source("scripts/data.R")

## Load model estimates
if (P$data_choice == "Vieira") {
  modfit <- readRDS(paste("mod_output/", P$data_choice, "_modfit.RDS", sep=""))
  
} else if (P$data_choice == "Mangiavacchi") {
  modfit <- readRDS(paste("mod_output/", P$data_choice, "_modfit.RDS", sep=""))
}

## Create lists to store plots (if not already in existence)
if(!exists("pv")) pv <- vector("list", length=2)
if(!exists("pm")) pm <- vector("list", length=3)


## (1) Non-stage-specific IgG
# Set up x-axis scale
if (P$data_choice == "Vieira") {
  x_axis   <- seq(10, 60, 10)
  x_limits <- c(10, 60)
  
} else if (P$data_choice == "Mangiavacchi") {
  x_axis   <- seq(0, 25, 5)
  x_limits <- c(0, 28.5)
}

# Set up colours, line types & point shapes
cols       <- c("mod_mean"="#000000", "mod_ribbon"="#999999", "dat" = "#999999")
line_types <- c("mod_mean"=1, "dat"=1)
shape      <- data.frame(c("mod_mean" = "a", "dat" = "b"))

# Make plot
p1 <- ggplot(data=modfit, aes(x=age, y=ig_mean)) +
  geom_point(aes(y=ig_mean, colour="mod_mean", fill="mod_mean", shape="b"), size=2) +
  geom_line(aes(y=ig_mean, group=1, colour="mod_mean", linetype="mod_mean"), size=0.5) +
  geom_ribbon(aes(ymin=ig_lower, ymax=ig_upper, fill="mod_ribbon"), alpha=0.30) +
  geom_point(data=df, aes(y=ig_mean, colour="dat", fill="dat", shape="a"), size=2) + 
  geom_errorbar(data=df, aes(ymin=ig_lower, ymax=ig_upper, group=2, colour="dat"), width=.5) +
  scale_colour_manual(name="", values=cols, guide=guide_legend(override.aes=aes(fill=NA))) + 
  scale_fill_manual(name="", values=cols, guide="none") + 
  scale_linetype_manual(name="", values=line_types) +
  labs(x="Age (years)", y="Seroprevalence IgG") + 
  scale_y_continuous(breaks=seq(0, 1, 0.20), limits = c(0,1), expand=c(0, 0)) + 
  scale_x_continuous(breaks=x_axis, limits=c(x_limits[1], x_limits[2]), expand=c(0, 0)) + 
  theme_light(base_size=9, base_line_size=0.5) +
  theme(legend.position='none') + 
  theme(axis.title=element_text(family="sans", size=9))

# Store plot in list
if (P$data_choice == "Vieira") {
  pv[[1]] <- p1
} else if (P$data_choice == "Mangiavacchi") {
  pm[[1]] <- p1
}


## (2) TgERP SD
p2 <- ggplot(data=modfit, aes(x=age, y=tgsd_mean)) +
  geom_line(aes(group=1, colour="mod_mean", linetype="mod_mean"), size=0.5) +
  geom_point(aes(colour="mod_mean", fill="mod_mean", shape="b"), size=2) +
  geom_ribbon(aes(ymin=tgsd_lower , ymax=tgsd_upper, fill="mod_ribbon"), alpha=0.20) +
  geom_point(data=df, aes(y=tgsd_mean, colour="dat", fill="dat", shape ="a"), size=2) + 
  geom_errorbar(data=df, aes(ymin=tgsd_lower, ymax=tgsd_upper, group=2, colour="dat"), width=.5) +
  scale_colour_manual(name="", values=cols, guide=guide_legend(override.aes=aes(fill=NA))) +
  scale_fill_manual(name="", values=cols, guide="none") +
  scale_linetype_manual(name="", values=line_types) +
  labs(x="Age (years)", y="Seroprevalence anti-TgERP IgG (SD)") + 
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1), expand=c(0, 0)) + 
  scale_x_continuous(breaks=x_axis, limits=c(x_limits[1], x_limits[2]), expand=c(0, 0)) + 
  theme_light(base_size=9, base_line_size=0.5) +
  theme(legend.position='none') + 
  theme(axis.title=element_text(family="sans", size=9))

# Store plot in list
if (P$data_choice == "Vieira") {
  pv[[2]] <- p2
} else if (P$data_choice == "Mangiavacchi") {
  pm[[2]] <- p2
}


## (3) TgERP ROC curve
if (P$data_choice == "Mangiavacchi") {
  p3 <- ggplot(data=modfit, aes(x=age, y=tgroc_mean)) +
    geom_line(aes(group=1, colour="mod_mean", linetype="mod_mean"), size=0.5) +
    geom_point(aes(colour="mod_mean", fill="mod_mean", shape="b"), size=2) +
    geom_ribbon(aes(ymin=tgroc_lower, ymax=tgroc_upper, fill="mod_ribbon"), alpha=0.20) +
    geom_point(data=df, aes(y=tgroc_mean, colour="dat", fill="dat", shape="a"), size=2) + 
    geom_errorbar(data=df, aes(ymin=tgroc_lower, ymax=tgroc_upper, group=2, colour="dat"), width=.5) +
    scale_colour_manual(name="", values=cols, guide=guide_legend(override.aes=aes(fill=NA))) +
    scale_fill_manual(name="", values=cols, guide="none") +
    scale_linetype_manual(name="", values=line_types) +
    labs(x="Age (years)", y="Seroprevalence anti-TgERP IgG (ROC)") + 
    scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1), expand=c(0, 0)) + 
    scale_x_continuous(breaks=x_axis, limits=c(x_limits[1], x_limits[2]), expand=c(0, 0)) + 
    theme_light(base_size=9, base_line_size=0.5) + 
    theme(legend.position='none') + 
    theme(axis.title=element_text(family="sans", size=9))
  
  # Store plot in list
  pm[[3]] <- p3
}


## Save multipanel plots of model-estimated seroprevalence vs. data
# Vieira et al. (2015)
ggarrange(pv[[1]], pv[[2]], nrow=2, font.label=list(size=11), labels=c("A", "B"))

ggsave(filename="plots/vieira_seroprev_multipanel.png",
       height=6, width=6, units="in", dpi=600)

# Mangiavacchi et al. (2016)
ggarrange(pm[[1]], ggarrange(pm[[3]], pm[[2]], ncol=2, font.label=list(size=11), labels=c("B", "C")), 
          nrow=2, labels="A", font.label=list(size=11))

ggsave(filename="plots/mangiavacchi_seroprev_multipanel.png",
       height=6, width=6, units="in", dpi=600)


## Plot estimates of anti-TgERP IgG duration
# Vieira et al. (2015)
tg_duration <- readRDS("mod_output/Vieira_sens1.RDS")
tg_duration <- tg_duration[c("tgerp", "tgerp_lower", "tgerp_upper")]  #extract cols needed 

# Make plot
p4 <- ggplot(data=tg_duration, aes(x=seq(1, 9), y=tgerp)) +
  geom_point(aes(colour="mod_mean"), size=2, group=2) +
  geom_errorbar(aes(ymin=tgerp_lower, ymax=tgerp_upper, group=2, colour="mod_mean"), width=.2) + 
  scale_colour_manual(name="",values=cols, guide = guide_legend(override.aes=aes(fill=NA))) + 
  labs(x="", y="Anti-TgERP IgG duration (years)") +
  scale_y_continuous(breaks=seq(0, 170, 20), limits=c(0,170), expand=c(0, 0)) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                     labels = c("0.80,\n0.85", "0.80,\n0.90", "0.80,\n0.95", 
                                "0.85,\n0.85", "0.85,\n0.90", "0.85,\n0.95", 
                                "0.90,\n0.85", "0.90,\n0.90", "0.90,\n0.95")) + 
  theme_light(base_size=9, base_line_size=0.5) + 
  theme(legend.position='none') + 
  theme(axis.title=element_text(family="sans", size=9))

# Mangiavacchi et al. (2016)
tg_duration <- readRDS("mod_output/Mangiavacchi_sens1.RDS")
tg_duration <- tg_duration[c("tgerp", "tgerp_lower", "tgerp_upper")]  #extract cols needed 

# Make plot
p5 <- ggplot(data=tg_duration, aes(x=seq(1, 9), y=tgerp)) +
  geom_point(aes(colour="mod_mean"), size=2, group=2) +
  geom_errorbar(aes(ymin=tgerp_lower, ymax=tgerp_upper, group=2, colour="mod_mean"), width=.2) + 
  scale_colour_manual(name="",values=cols, guide = guide_legend(override.aes=aes(fill=NA))) + 
  labs(x="", y="Anti-TgERP IgG duration (years)") +
  scale_y_continuous(breaks=seq(0, 170, 20), limits=c(0,170), expand=c(0, 0)) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                     labels = c("0.80,\n0.85", "0.80,\n0.90", "0.80,\n0.95", 
                                "0.85,\n0.85", "0.85,\n0.90", "0.85,\n0.95", 
                                "0.90,\n0.85", "0.90,\n0.90", "0.90,\n0.95")) + 
  theme_light(base_size=9, base_line_size=0.5) + 
  theme(legend.position='none') + 
  theme(axis.title=element_text(family="sans", size=9))

## Save multipanel plot of anti-TgERP IgG duration estimates
ggarrange(p4, p5, nrow=2, font.label=list(size=11), labels=c("A", "B"))

ggsave(filename="plots/TgERP_duration_multipanel.png",
       height=6, width=6, units="in", dpi=600)

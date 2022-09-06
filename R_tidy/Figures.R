#
## Figures
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#


# Load packages ----

packages <- c("ggplot2", "ggpubr")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load results ----

load("Results/single_plots.rda")


## Figure 5 ----

Fig5_Reg <- ggarrange(cl_IS_TL + rremove("ylab"), cl_IS_DEG + rremove("ylab"), 
                      cl_IS_TS + rremove("ylab"), cl_IS_HAB + rremove("ylab"), 
                      common.legend = TRUE, legend = "bottom",
                      labels = c("A", "B", "C", "D"),
                      ncol=2, nrow=2)
Fig5_Reg <- annotate_figure(Fig5_Reg, left = text_grob("log(mean Interaction Strength)", rot = 90, 
                                                       vjust = 1, size = 18))
Fig5_Reg

# ggsave(filename = "Manuscript/Fig.5_Reg.png", plot = Fig5_Reg, 
#        width = 10, units = "in", dpi = 600, bg = "white")



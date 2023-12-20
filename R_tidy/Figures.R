#
## Figures
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022-2023
#


# Load packages -----------------------------------------------------------
packages <- c("ggplot2", "ggpubr", "multiweb", "igraph", "factoextra", "cluster", "scales")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load results ------------------------------------------------------------

load("Results/interaction_estimation_sim.rda")
load("Results/net_&_spp_prop_sim.rda")
load("Results/single_plots_sep23.rda")
load("Results/QSS_summary_oct31.rda")


# Figures -----------------------------------------------------------------
## Figure 1 ----
# Map of the Weddell Sea and Dronning Maud Land sector highlighting the high Antarctic shelf as a dashed-line contour.
# Modified from www.soos.aq.

## Figure 2 ----
# Scheme of a network showing the weighted and unweighted properties we used to characterize the species 
# of the Weddell Sea food web

## Figure 3 ----
# Frequency distribution of interaction strengths for the Weddell Sea food web (n = 490).
# IS distribution
# options(scipen = 999)  # avoid scientific notation
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

Fig3_IntDist <- ggplot(wedd_int_pd_summary, aes(IS_med)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  scale_y_log10() +
  scale_x_continuous(label = scientific_10) +
  labs(x = "Interaction Strength (median)", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))
Fig3_IntDist

ggsave(filename = "Manuscript/Fig3_IntDist_sim.png", plot = Fig3_IntDist,
       width = 10, units = "in", dpi = 600, bg = "white")


## Figure 4 ----
# Relationships between weighted (interaction strength IS) and unweighted properties including habitat.
# Linear regressions are shown between log(mean interaction strength) and trophic level (A), degree (B) 
# and trophic similarity (C) for each group of species resulted from the clustering analysis: 'High IS' and 'Low IS'.
# All regression slopes are significant and statistically different among groups (p < 0.01).

Fig4_LinReg <- ggarrange(cl_IS_TL + rremove("ylab"), cl_IS_DEG + rremove("ylab"), 
                      cl_IS_TS + rremove("ylab"), cl_IS_HAB + rremove("ylab"), 
                      common.legend = TRUE, legend = "bottom",
                      labels = c("A", "B", "C", "D"),
                      ncol=2, nrow=2)
Fig4_LinReg <- annotate_figure(Fig4_LinReg, left = text_grob("log(median Interaction Strength)", rot = 90, 
                                                       vjust = 1, size = 18))
Fig4_LinReg

ggsave(filename = "Manuscript/Fig4_LinReg_sim.png", plot = Fig4_LinReg,
       width = 10, units = "in", dpi = 600, bg = "white")


## Figure 5 ----
# Quasi-Sign Stability (QSS) difference between the whole Weddell Sea food web (n = 490) and the food web 
# without one species (n = 489) for weighted (interaction strength) and unweighted species properties, and habitat.
# Each point represents a species colored by group ('High IS' or 'Low IS'). Shape indicates the impact on the QSS; 
# if significant the extinction of that species altered the stability of the food web.

legend <- ggplot(all_dif, aes(x = meanTrophicSimil, y = difQSSrelat)) +
  geom_point(aes(color = ifelse(coding>0, "> 0.55", "< 0.55"))) + 
  scale_color_viridis_d(direction=-1) +
  # ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS))+
  # geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  # scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  lims(x = c(0,0), y = c(0,0)) +
  theme_void() +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  10),
        legend.title = element_text(size = 12, face = "bold")) +
  guides(color = guide_legend(title = "Stability impact"))

Fig.5_QSSDif <- ggarrange(IS_QSS + rremove("ylab"), TL_QSS + rremove("ylab"), 
                      DEG_QSS + rremove("ylab"), TS_QSS + rremove("ylab"), 
                      HAB_QSS + rremove("ylab"), legend,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol=3, nrow=2)
Fig.5_QSSDif <- annotate_figure(Fig.5_QSSDif, left = text_grob("Stability difference", rot = 90, 
                                                       vjust = 1, size = 18))
Fig.5_QSSDif

ggsave(filename = "Manuscript/Fig.5_QSSDif.png", plot = Fig.5_QSSDif,
       width = 10, units = "in", dpi = 600, bg = "white")


# Appendix ----------------------------------------------------------------
## App 1 ----
# Graphic representation of the Weddell Sea food web. Species (nodes) are arranged vertically and colored by trophic level.
# The diameter of the node indicates the total number of interactions. Predator-prey interactions are represented by the arrows,
# from prey to predator.

# Food web
App1_FWplot <- plot_troph_level(g, vertexSizeFactor = V(g)$Deg*0.04, edge.width = 0.25,
                                ylab = "Trophic level")

# ggsave(filename = "Manuscript/App1_FWplot.png", plot = App1_FWplot,
#        width = 10, units = "in", dpi = 600, bg = "white")
# Error in UseMethod("grid.draw") : 
# no applicable method for 'grid.draw' applied to an object of class "c('double', 'numeric')"


# Supplementary Material --------------------------------------------------

## Figure S1 ----
# Frequency distribution of interquartile range for the estimated interaction strengths (n=16041) of the Weddell Sea food web
IQR <- ggplot(wedd_int_pd_summary, aes(IS_iqr)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  scale_y_log10() +
  labs(x = "Interquartil of IS estimation", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

ggsave(filename = "Manuscript/Supp1_IQR.png", plot = IQR,
       width = 10, units = "in", dpi = 600, bg = "white")

## Figure S2 ----
# Distribution of relative stability differences (between the whole network and the network minus one species) when the species in question are removed from the Weddell Sea food web. Central tendencies are shown: median in brown dash, mode in blue longdash.
QSS_distr
ggsave(filename = "Manuscript/Supp2_QSS_distr.png", 
       width = 10, units = "in", dpi = 600, bg = "white")

## Table S1 ----
# Weighted (interaction strength) and unweighted properties of the trophic species of Weddell Sea food web. Ordered by decreasing median interaction strength. median IS = median interaction strength, Q1 IS = First quartil of the IS distribution, Q3 IS = Third quartil of the IS distribution, TL = trophic level, TS = trophic similarity.
names(all_dif)
knitr::kable(all_dif %>% dplyr::filter(coding==1) %>% 
               dplyr::select(TrophicSpecies, IS_median, IS_Q1, IS_Q3, TotalDegree, TL, meanTrophicSimil, median_difQSS) %>% 
               dplyr::arrange(., desc(IS_median)))

## Table S2 ----
# Summary of maximum eigenvalue (QSS) distribution of differences before and after performing extinction simulations in the Weddell Sea food web. Ordered by decreasing proportion of positive differences. Prop dif QSS +  = Proportion of positive differences, Prop dif QSS - = Proportion of negative differences, median difQSS relat = median of relative QSS differences.
knitr::kable(all_dif %>% 
               dplyr::select(TrophicSpecies, prop_difQSSm_pos, prop_difQSSm_neg, median_difQSS) %>% 
               dplyr::arrange(., prop_difQSSm_neg))

#
## Figures
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#


# Load packages ----

packages <- c("ggplot2", "ggpubr", "multiweb", "igraph", "factoextra", "cluster")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load results ----

load("Results/interaction_estimation.rda")
load("Results/net_&_spp_prop.rda")
load("Results/single_plots.rda")
#load("Results/QSS_extinction_dif.rda")
#load("Results/QSS_summary.rda")


# Figures ----

## Figure 1 ----
# Map of the Weddell Sea and Dronning Maud Land sector highlighting the high Antarctic shelf as a dashed-line contour.
# Modified from www.soos.aq.

## Figure 2 ----
# Scheme of a network showing the weighted and unweighted properties we used to characterize the species 
# of the Weddell Sea food web


## Figure 3 ----
# Frequency distribution of interaction strengths for the Weddell Sea food web (n = 490).

# Interaction strength distribution
Fig3_IntDist <- ggplot(wedd_int_pd, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  scale_y_log10() +
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))
Fig3_IntDist

ggsave(filename = "Manuscript/Fig3_IntDist.png", plot = Fig3_IntDist,
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
Fig4_LinReg <- annotate_figure(Fig4_LinReg, left = text_grob("log(mean Interaction Strength)", rot = 90, 
                                                       vjust = 1, size = 18))
Fig4_LinReg

ggsave(filename = "Manuscript/Fig4_LinReg.png", plot = Fig4_LinReg,
       width = 10, units = "in", dpi = 600, bg = "white")


## Figure 5 ----
# Quasi-Sign Stability (QSS) difference between the whole Weddell Sea food web (n = 490) and the food web 
# without one species (n = 489) for weighted (interaction strength) and unweighted species properties, and habitat.
# Each point represents a species colored by group ('High IS' or 'Low IS'). Shape indicates the impact on the QSS; 
# if significant the extinction of that species altered the stability of the food web.

legend <- ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS))+
  geom_point(aes(color = cluster_mean, 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  lims(x = c(0,0), y = c(0,0)) +
  theme_void() +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  10),
        legend.title = element_text(size = 12, face = "bold")) +
  guides(color = guide_legend(title = "Group"), shape = guide_legend(title = "Stability impact"))

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


# Appendix ----
## App 1 ----
# Graphic representation of the Weddell Sea food web. Species (nodes) are arranged vertically and colored by trophic level.
# The diameter of the node indicates the total number of interactions. Predator-prey interactions are represented by the arrows,
# from prey to predator.

# Food web
App1_FWplot <- plot_troph_level(g, vertexSizeFactor = V(g)$Deg*0.04, edge.width = 0.25)

# ggsave(filename = "Manuscript/App1_FWplot.png", plot = App1_FWplot,
#        width = 10, units = "in", dpi = 600, bg = "white")
# Error in UseMethod("grid.draw") : 
# no applicable method for 'grid.draw' applied to an object of class "c('double', 'numeric')"


## App 2 ----
# A. Frequency distribution for the mean interaction strength of the species of Weddell Sea food web. 
# B. Visualization of the optimal number of clusters applying the Gap statistic.

a <- ggplot(cluster_data, aes(x = log(IS_mean), color = cluster_mean, fill= cluster_mean)) + 
  geom_density(alpha=0.3) + 
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_fill_manual(values = c("#541352FF", "#ffcf20FF")) +
  labs(x = "log(mean Interaction Strength)", y = "Frequency", color = "Cluster") +
  theme_bw() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
a <- a + guides(fill = "none")

# Determine and visualize the optimal number of clusters using total within sum of square
data.scaled <- scale(log(spp_all_prop$IS_mean))
fviz_nbclust(data.scaled, kmeans, method = "wss")
# Calculate gap statistic based on number of clusters
gap_stat <- clusGap(data.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
# Plot number of clusters vs. gap statistic
b <- fviz_gap_stat(gap_stat) +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
b

App2_Clus_meanIS <- ggarrange(a, b,
                         labels = c("A", "B"),
                         ncol=2)
App2_Clus_meanIS

ggsave(filename = "Manuscript/App2_Clus_meanIS.png", plot = App2_Clus_meanIS,
       width = 10, units = "in", dpi = 600, bg = "white")


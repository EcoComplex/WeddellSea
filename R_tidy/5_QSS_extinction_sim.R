#
## Calculate Quasi-sign stability 'QSS' and compare before and after extinction simulations
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first run 'R/1_calc_interaction_strength.R', 
# 'R/3_species_w&uw_prop.R' & 'R/4_quantile_regression.R'


# Load packages ----

packages <- c("igraph", "multiweb", "dplyr", "tictoc", "ggplot2")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

load("Results/net_&_spp_prop.rda")


# Run extinction simulations & estimate QSS ----

# List of species to delete
sp_list <- arrange(spp_all_prop, desc(IS_mean)) %>% 
  dplyr::select(TrophicSpecies, IS_mean)
sp_list <- sp_list$TrophicSpecies

# Be aware that this simulation might take days
# To reduce the time you can configure the number of simulations 'nsim'
nsim <- 1000
print(paste("QSS 1 sp extinction  - Nsim = ",nsim))
tic("QSS dif")
QSS_extinction_dif <- calc_QSS_extinction_dif(g, sp_list, ncores=48, nsim=nsim, istrength = TRUE)
toc()

QSS_extinction_dif <- as_tibble(QSS_extinction_dif)

## Save simulations ----

save(QSS_extinction_dif,
     file = "Results/QSS_extinction_dif.rda")


# QSS vs spp prop ----

# Load needed data
load("Results/QSS_extinction_dif.rda")

all_data <- QSS_extinction_dif %>% 
  rename(TrophicSpecies = Deleted) %>% 
  left_join(spp_all_prop)

# Species w/ QSS significant impact
QSS_sig <- all_data %>% 
  dplyr::filter(Ad_pvalue < 0.01) %>% 
  dplyr::select(TrophicSpecies, IS_mean, TL, TotalDegree, meanTrophicSimil, Habitat, difQSS, Ad_pvalue) %>% 
  arrange(., Ad_pvalue)


## By mean interaction strength ----
IS_QSS <- ggplot(all_data, aes(x = log(IS_mean), y = difQSS)) +
  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Stability impact", x = "log(mean Interaction Strength)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
IS_QSS

## By trophic level ----
TL_QSS <- ggplot(all_data, aes(x = TL, y = difQSS)) +
  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(shape = "Stability impact", x = "Trophic level", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
TL_QSS

## By degree ----
DEG_QSS <- ggplot(all_data, aes(x = TotalDegree, y = difQSS)) +
  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  scale_x_log10() +
  labs(color = "Stability impact", x = "Degree (log scale)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
DEG_QSS

## By trophic similarity ----
TS_QSS <- ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS)) +
  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Stability impact", x = "Trophic similarity", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
TS_QSS

## By habitat ----
HAB_QSS <- ggplot(all_data, aes(x = Habitat, y = difQSS)) +
  geom_violin(fill = "grey90") +
  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Stability impact", x = "Habitat", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 12))
HAB_QSS


# Save results ----

save(all_data, IS_QSS, TL_QSS, DEG_QSS, TS_QSS, HAB_QSS,
     file = "Results/QSS_summary_sep22.rda")

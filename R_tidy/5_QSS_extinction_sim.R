#
## Calculate Quasi-sign stability 'QSS' and compare before and after extinction simulations
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first run 'R/1_calc_interaction_strength.R', 
# 'R/3_species_w&uw_prop.R' & 'R/4_quantile_regression.R'



# Load packages -----------------------------------------------------------


packages <- c("igraph", "dplyr", "tictoc", "ggplot2", "devtools")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)

if( require(multiweb) == FALSE) {
  devtools::install_github("lsaravia/multiweb")
  require(multiweb)
}

# Load data ---------------------------------------------------------------
load("Results/net_&_spp_prop_sim.rda")

# Extinction simulations & QSS --------------------------------------------
# List of species to delete
sp_list <- arrange(spp_all_prop, desc(IS_mean)) %>% 
  dplyr::select(TrophicSpecies, IS_mean)
sp_list <- sp_list$TrophicSpecies

# Be aware that this simulation might take days
# To reduce the time you may configure the number of simulations 'nsim'
nsim <- 1000
print(paste("QSS 1 sp extinction SIM 3 - Nsim = ", nsim))
tic("QSS dif")
QSS_extinction_dif <- calc_QSS_extinction_dif(g, sp_list, ncores=48, nsim=nsim, istrength = TRUE)
toc()

QSS_extinction_dif <- as_tibble(QSS_extinction_dif)

## Save simulations ----
# save(QSS_extinction_dif,
#      file = "Results/QSS_extinction_dif.rda")


# QSS vs spp prop ----
# Load previous results 
load("Results/QSS_summary_oct30.rda")

# bind with new simulation
#
all_data_new_sim <- QSS_extinction_dif %>% 
  rename(TrophicSpecies = Deleted) %>% 
  left_join(spp_all_prop) %>% mutate(sim=3)       # Simulation number 3

all_data <- bind_rows(all_data, all_data_new_sim)

# Check the species that appear in all simulations and filter by difQSS 1% 
QSS_sig <- all_data %>%
  group_by(TrophicSpecies) %>%
  filter(n() > 1 & all(difQSS > 0) | all(difQSS < 0)) %>% slice_max(order_by = abs(difQSS), n = 1) %>%
  ungroup() %>% mutate(difQSSrelat = difQSS/QSS_all) %>% filter(abs(difQSSrelat) > 0.01) %>%
  dplyr::select(TrophicSpecies, IS_mean, TL, TotalDegree, meanTrophicSimil, Habitat, difQSS,difQSSrelat, Ad_pvalue) %>% 
  arrange(., difQSSrelat)  


# Species w/ QSS significant impact
# QSS_sig <- all_data %>% 
#   dplyr::filter(Ad_pvalue < 0.05) %>% mutate(difQSSrelat = difQSS/QSS_all) %>%
#   dplyr::select(TrophicSpecies, IS_mean, TL, TotalDegree, meanTrophicSimil, Habitat, difQSS,difQSSrelat, Ad_pvalue) %>% 
#   arrange(., Ad_pvalue)


## By mean interaction strength ----
IS_QSS <- ggplot(all_data, aes(x = log(IS_mean), y = difQSS,text=TrophicSpecies)) +
  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) + scale_color_viridis_d(direction=-1) +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
#  geom_point(aes(color = Ad_pvalue)) + scale_color_viridis_c() +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Stability impact", x = "log(mean Interaction Strength)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
IS_QSS
ggplotly(IS_QSS, tooltip=c("x", "y", "text"))

## By trophic level ----
TL_QSS <- ggplot(all_data, aes(x = TL, y = difQSS,text=TrophicSpecies)) +
#  geom_point(aes(color = ifelse(Ad_pvalue < 0.05, "Significant", "Non-significant"))) +
  geom_point(aes(color = Ad_pvalue)) + scale_color_viridis_c() +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  
  labs(shape = "Stability impact", x = "Trophic level", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
TL_QSS
require(plotly)
ggplotly(TL_QSS, tooltip=c("x", "y", "text"))

## By degree ----
DEG_QSS <- ggplot(all_data, aes(x = TotalDegree, y = difQSS,text=TrophicSpecies)) +
#  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  geom_point(aes(color = Ad_pvalue)) + scale_color_viridis_c() +
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
ggplotly(DEG_QSS, tooltip=c("x", "y", "text"))

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
  geom_jitter(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant")),width = 0.05,alpha=0.5) +
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


# Save results ------------------------------------------------------------

#save(all_data, IS_QSS, TL_QSS, DEG_QSS, TS_QSS, HAB_QSS,
#     file = "Results/QSS_summary_sep22.rda")

save(all_data, IS_QSS, TL_QSS, DEG_QSS, TS_QSS, HAB_QSS,
          file = "Results/QSS_summary_oct31.rda")

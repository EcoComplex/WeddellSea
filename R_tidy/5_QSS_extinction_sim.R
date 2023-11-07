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
nsim <- 10
print(paste("QSS 1 sp extinction NEW version - Nsim = ", nsim))
tic("QSS dif")
QSS_extinction_dif <- calc_QSS_extinction_dif(g, sp_list, ncores=48, nsim=nsim, istrength = TRUE)
toc()

QSS_extinction_dif <- as_tibble(QSS_extinction_dif)

## Save simulations ----
save(QSS_extinction_dif,
      file = "Results/QSS_extinction_dif.rda")


# QSS vs spp prop ----
# Load previous results 
#
# Need to change this --- no exec in the server
#
if(FALSE){

load("Results/QSS_extinction_dif.rda")

# * We need to find the important species
  
# * Then summarize the medians
  
# Calculate also the diffQSS from the median QSS_all
#
all_data_new <- QSS_extinction_dif %>% 
  rename(TrophicSpecies = Deleted) %>% group_by(TrophicSpecies) %>% 
  mutate(difQSSm = median(QSS_all) - QSS_ext, difQSSrelat = difQSSm/QSS_all) # %>% left_join(spp_all_prop) 

# Calculate the proportion of >0 or <0
#
props <- all_data_new %>%
  group_by(TrophicSpecies) %>%
  summarize(
    prop_difQSS_pos = mean(difQSS > 0),
    prop_difQSS_neg = mean(difQSS < 0),
    prop_difQSSm_pos = mean(difQSSm > 0), 
    prop_difQSSm_neg = mean(difQSSm < 0),
  ) 

all_data_new <- all_data_new %>%
  left_join(props)

# Define a function to estimate the mode from a density object
estimate_mode <- function(x) {
  dens <- density(x)
  x_value <- dens$x[which.max(dens$y)]
  return(x_value)
}

# Calculate the mode for each TrophicSpecies
mode_values <- all_data_new %>%
  group_by(TrophicSpecies) %>%
  summarize(
    mode_difQSS = estimate_mode(difQSSrelat),
    median_difQSS = median(difQSSrelat)
  )

# Join mode values back to the main data frame
all_data_new <- all_data_new %>%
  left_join(mode_values, by = "TrophicSpecies")


# Plot with the proportion
#
QSS_distr <- all_data_new %>%
  filter(TrophicSpecies %in% c("Orcinus orca","Hydrurga leptonyx","Euphausia superba","Balaenoptera acutorostrata")) %>% 
  ggplot(aes(difQSSrelat,fill=TrophicSpecies)) + geom_density( alpha=0.3 ) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  facet_wrap(~TrophicSpecies, ncol = 2) + scale_fill_viridis_d(guide=FALSE) +
  geom_text(aes(x = 0, y = 0.5, label = paste0( prop_difQSSm_neg*100,"%")), 
            position = position_nudge(x = -0.2)) +
  geom_text(aes(x = 0, y = 0.5, label = paste0( prop_difQSSm_pos*100,"%")),
            position = position_nudge(x = 0.2)) + 
  geom_vline( aes(xintercept = mode_difQSS), color = "blue", linetype = "longdash") +
  geom_vline( aes(xintercept = median_difQSS), color = "brown", linetype = "dashed") +
  labs(x = "Stability difference %", y = "Density") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))

QSS_distr

# Check the species that appear in all simulations and filter by difQSS 1% 
QSS_sig <- props %>%
  filter(prop_difQSSm_pos > 0.55 | prop_difQSSm_neg > 0.55) 

all_data_TS <- all_data_new %>% group_by(TrophicSpecies) %>% 
  summarize(
    prop_difQSS_pos = mean(difQSS > 0),
    prop_difQSS_neg = mean(difQSS < 0),
    prop_difQSSm_pos = mean(difQSSm > 0), 
    prop_difQSSm_neg = mean(difQSSm < 0),
    difQSS=mean(difQSS),
    mode_difQSS = estimate_mode(difQSSrelat),
    median_difQSS = median(difQSSrelat),
    difQSSrelat=mean(difQSSrelat)
    ) %>% left_join(spp_all_prop) 


# Species w/ QSS significant impact
# 
all_dif <- all_data_TS %>%
  left_join(QSS_sig %>% select(TrophicSpecies,prop_difQSS_pos) ,by="TrophicSpecies") %>% 
  mutate(coding=ifelse(!is.na(prop_difQSS_pos.y),ifelse(difQSSrelat>0,1,1),0))

## By median interaction strength ----
IS_QSS <- ggplot(all_dif, aes(x = log(IS_median), y = difQSSrelat,text=TrophicSpecies)) +
  geom_point(aes(color = coding),alpha=0.6) + scale_color_viridis_c(direction=-1) +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
#  geom_point(aes(color = Ad_pvalue)) + scale_color_viridis_c() +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Stability impact", x = "log(median Interaction Strength)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
IS_QSS
#require(plotly)
#ggplotly(IS_QSS, tooltip=c("x", "y", "text"))

## By trophic level ----
TL_QSS <- ggplot(all_dif, aes(x = TL, y = difQSSrelat,text=TrophicSpecies)) +
  geom_point(aes(color = coding),alpha=0.6) + scale_color_viridis_c(direction=-1) +
#  geom_point(aes(color = ifelse(Ad_pvalue < 0.05, "Significant", "Non-significant"))) +
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
#ggplotly(TL_QSS, tooltip=c("x", "y", "text"))

## By degree ----
DEG_QSS <- ggplot(all_dif, aes(x = TotalDegree, y = difQSSrelat,text=TrophicSpecies)) +
#  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  geom_point(aes(color = coding),alpha=0.6) + scale_color_viridis_c(direction=-1) +
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
#ggplotly(DEG_QSS, tooltip=c("x", "y", "text"))

## By trophic similarity ----
TS_QSS <- ggplot(all_dif, aes(x = meanTrophicSimil, y = difQSSrelat)) +
  geom_point(aes(color = coding),alpha=0.6) + scale_color_viridis_c(direction=-1) +
#  geom_point(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
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
(HAB_QSS <- ggplot(all_dif, aes(x = Habitat, y = difQSSrelat)) +
  geom_violin(fill = "grey90") +
  geom_jitter(aes(color = coding),width=0.08,alpha=0.6) + 
   geom_point(data=all_dif %>% filter(Habitat=="Benthic" & coding==1),alpha=0.6) +
   scale_color_viridis_c(direction=-1) +
#  geom_jitter(aes(color = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant")),width = 0.05,alpha=0.5) +
#  scale_color_manual(values = c("black", "red"), labels = c("Non-significant", "Significant")) +
  #scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Stability impact", x = "Habitat", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 12))
)


# Save results ------------------------------------------------------------

#save(all_data, IS_QSS, TL_QSS, DEG_QSS, TS_QSS, HAB_QSS,
#     file = "Results/QSS_summary_sep22.rda")

save(all_data_new, all_dif, IS_QSS, TL_QSS, DEG_QSS, TS_QSS, HAB_QSS, QSS_distr,
          file = "Results/QSS_summary_oct31.rda")
}


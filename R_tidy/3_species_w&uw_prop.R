#
## Calculate weighted (interaction strength) & unweighted species properties
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first run 'R/1_calc_interaction_strength.R'


# Load packages ----

packages <- c("tidyverse", "igraph", "NetIndices", "cheddar", "multiweb")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Weighted properties ----

## Load data ----
load("Results/interaction_estimation_sim.rda")


## Interaction strength per species ----
### All IS ----
# Need to run this before:
wedd_int_pd <- multiweb::calc_interaction_intensity(wedd_df_pd, res.mass.mean.kg., res_den, con.mass.mean.kg., interaction.dimensionality, nsims=1000)
# Interaction strength for Consumers (incoming)
IS_con_all <- wedd_int_pd %>% 
  dplyr::select(con.taxonomy, qRC) %>% 
  dplyr::rename(TrophicSpecies = con.taxonomy)
# Interaction strength for Resources (outcoming)
IS_res_all <- wedd_int_pd %>% 
  dplyr::select(res.taxonomy, qRC) %>% 
  dplyr::rename(TrophicSpecies = res.taxonomy)
# Bind interaction strengths for each species
spp_w_prop <- bind_rows(IS_con_all, IS_res_all) %>% 
  group_by(TrophicSpecies) %>% 
  summarise(ci = list(enframe(Hmisc::smean.cl.normal(qRC)))) %>% 
  unnest %>% 
  spread(name, value) %>% 
  rename(IS_low_ci=Lower, IS_high_ci=Upper, IS_mean=Mean)


### Mean IS ----
# Select interaction strength for Consumers (incoming)
IS_con_m <- wedd_int_pd_summary %>% 
  dplyr::select(con.taxonomy, IS_mean) %>% 
  dplyr::rename(TrophicSpecies = con.taxonomy)
# Select interaction strength for Resources (outcoming)
IS_res_m <- wedd_int_pd_summary %>% 
  dplyr::select(res.taxonomy, IS_mean) %>% 
  dplyr::rename(TrophicSpecies = res.taxonomy)
# Bind interaction strengths for each species
# Calculate mean, Q1 & Q3
IS_mean <- bind_rows(IS_con_m, IS_res_m)
IS_all_m <- bind_rows(IS_con_m, IS_res_m) %>% 
  dplyr::group_by(TrophicSpecies) %>% 
  dplyr::summarize(IS_mean = mean(IS_mean),
                   IS_median = median(IS_mean),
                   IS_Q1 = quantile(IS_mean, 0.25),
                   IS_Q3 = quantile(IS_mean, 0.75),
                   IS_max = max(IS_mean),
            Check_NumbInt = n())

# Convert interaction list to an igraph with weights
g_mean <- graph_from_data_frame(wedd_int_pd_summary %>% 
                             dplyr::select(res.taxonomy, con.taxonomy, IS_mean) %>% 
                             dplyr::rename(weight=IS_mean), directed = TRUE)
E(g_mean)$weight
# Add sum of in & out interaction strengths per species
# V(g)$IS_sum_in <- strength(g, mode = "in")
# V(g)$IS_sum_out <- strength(g, mode = "out")
# V(g)$IS_sum_tot <- strength(g, mode = "total")
# vertex_attr_names(g)


# Unweighted properties ----
# Trophic level
adj_mat <- as_adjacency_matrix(g_mean, sparse = TRUE)
tl <- round(TrophInd(as.matrix(adj_mat)), digits = 3)
V(g_mean)$TL <- tl$TL
V(g_mean)$Omn <- tl$OI
vertex_attr_names(g_mean)
# Degree
V(g_mean)$Deg <- degree(g_mean, mode = "total")
V(g_mean)$InDeg <- degree(g_mean, mode = "in")
V(g_mean)$OutDeg <- degree(g_mean, mode = "out")
vertex_attr_names(g_mean)
# Trophic similarity
source("R/igraph_cheddar.R")  # load function to convert igraph to cheddar object
igraph_to_cheddar(g_mean)
cc <- LoadCommunity("Community")
ts <- TrophicSimilarity(cc)
mts <- tibble(TrophicSpecies = rownames(ts), meanTrophicSimil = colMeans(ts))  # data frame


# All spp properties ----
# Weighted
# spp_IS_sum_in <- data.frame(V(g)$name, V(g)$IS_sum_in)
# colnames(spp_IS_sum_in) <- c("TrophicSpecies", "IS_sum_in")
# spp_IS_sum_out <- data.frame(V(g)$name, V(g)$IS_sum_out)
# colnames(spp_IS_sum_out) <- c("TrophicSpecies", "IS_sum_out")
# spp_IS_sum_tot <- data.frame(V(g)$name, V(g)$IS_sum_tot)
# colnames(spp_IS_sum_tot) <- c("TrophicSpecies", "IS_sum_tot")
# spp_w_prop <- total_int %>% 
#   left_join(spp_IS_sum_in) %>% 
#   left_join(spp_IS_sum_out) %>% 
#   left_join(spp_IS_sum_tot)

# Unweighted
spp_name <- as.data.frame(V(g_mean)$name)
spp_totdeg <- as.data.frame(V(g_mean)$Deg)
spp_indeg <- as.data.frame(V(g_mean)$InDeg)
spp_outdeg <- as.data.frame(V(g_mean)$OutDeg)
spp_tl <- as.data.frame(V(g_mean)$TL)
spp_omn <- as.data.frame(V(g_mean)$Omn)
spp_tropsimil <- mts[,2]
spp_uw_prop <- bind_cols(spp_name, spp_totdeg, spp_indeg, spp_outdeg, spp_tl,
                         spp_omn, spp_tropsimil)
colnames(spp_uw_prop) <- c("TrophicSpecies", "TotalDegree", "InDegree", "OutDegree",
                        "TL", "Omn", "meanTrophicSimil")

# Bind weighted & unweighted properties
spp_all_prop <- spp_uw_prop %>% 
  left_join(spp_w_prop)

# Load & include Habitat
hab_data <- read.csv(file = "Data/WeddellSeaHabitat.csv")
hab_data_inc <- hab_data %>% 
  dplyr::select(species, Environment) %>% 
  dplyr::rename(TrophicSpecies = species, Habitat = Environment)
hab_data_comp <- hab_data_inc %>% 
  add_row(TrophicSpecies = "Phytodetritus", Habitat = "Benthic") %>%            # missing species
  add_row(TrophicSpecies = "Sediment", Habitat = "Benthic") %>%                 # missing species
  add_row(TrophicSpecies = "Iphimediella  cyclogena", Habitat = "Benthic") %>%  # missing species
  mutate(Habitat = case_when(Habitat == "Bathydemersal" ~ "Demersal", TRUE ~ Habitat))  # replace 'Bathydemersal' with 'Demersal'

spp_all_prop <- spp_all_prop %>% 
  left_join(hab_data_comp)


# Save results ----

save(g_mean, IS_all_m, spp_w_prop, spp_uw_prop, spp_all_prop,
     file = "Results/net_&_spp_prop_sim.rda")

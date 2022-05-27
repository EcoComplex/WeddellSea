## Comparing unweighted and weighted results
# using 'multiweb' package & Susanne Kortsch's function 'fluxind'


# Load packages ----

require(igraph)
require(multiweb)
require(tidyverse)
# Load function to estimate interaction intensity & 'fluxind'
source("network_fun.r")

# Load data ----

# Original data (from GATEWAy)
ga <- readr::read_csv("../Data/283_2_FoodWebDataBase_2018_12_10.csv",col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
names(ga)
wedd_df <-  ga %>%filter(grepl('Weddell',foodweb.name))

# Tidy data
wedd_df <- wedd_df %>%
  rename(con_mass_mean=con.mass.mean.g., res_mass_mean=res.mass.mean.g., con_taxonomy=con.taxonomy,
         res_taxonomy=res.taxonomy, interaction_dim=interaction.dimensionality) %>%
  dplyr::select(con_taxonomy, res_taxonomy, con_mass_mean,res_mass_mean, interaction_dim)

# Complete interaction dimensionality data
wedd_df <- read_delim("../Data/Wedd_mass_complete.dat", delim=" ")

# Unweighted food web
g_uw <- graph_from_data_frame(wedd_df %>% select(res_taxonomy, con_taxonomy),
                              directed=TRUE)


# Calculate interaction intensity ----
# Using 'interaction_intensity' function from 'network_fun.r'

wedd_int <- interaction_intensity(wedd_df, res_mass_mean, con_mass_mean,
                                  interaction_dim)

# Explore distribution of interaction intensity
ggplot(wedd_int, aes(qRC)) +
  geom_histogram(bins=50, color="darkblue", fill="white") +
  theme_bw() +
  scale_y_log10()

# Weighted food web
g_w <- graph_from_data_frame(wedd_int %>%
                             select(res_taxonomy, con_taxonomy, qRC) %>%
                             rename(weight=qRC), directed=TRUE)
E(g_w)$weight


# Calculate properties ----

# Multiweb package
prop_mw <- calc_topological_indices(g_uw)

# Susanne Kortsch functions
adj_uw <- as_adjacency_matrix(g_uw, sparse = TRUE)
adj_w <- as_adjacency_matrix(g_w, attr = "weight", sparse = TRUE)

# How can the unweighted results be different?
prop_sk_uw <- fluxind(as.matrix(adj_uw))
prop_sk_w <- fluxind(as.matrix(adj_w))

# # Trophic level
# tl_uw <- TrophInd(as.matrix(adj_uw))
# tl_w <- TrophInd(as.matrix(adj_w))
# 
# indegree <- degree(g_uw, mode = "in")
# V(g_uw)$indeg <- indegree
# vertex.attributes(g_uw)
# 
# TroLev <- function(fw)
# {
#   fw <- t(fw)
#   nn <- rowSums(fw); nn[nn==0] <- 1
#   ww <- diag(1/nn)
#   L1 <- ww %*% fw
#   L2 <- L1 - diag(rep(1,length(nn)))
#   b <- -1*rep(1,length(nn))
#   Tro.lev <- solve(L2) %*% b
#   
#   return(Tro.lev)
# }
# TroLev(as.matrix(adj_w))
# 
# # Generality & Vulnerability
# gen.fun <- function(g){
#   V(g)$indeg <- degree(g, mode = "in")
#   data_indeg <- as.data.frame(V(g)$indeg) %>% 
#     filter(V(g)$indeg != 0) %>% 
#     rename(n = "V(g)$indeg")
#   G_avg <- round(mean(data_indeg$n), 2)
#   G_med <- round(median(data_indeg$n), 2)
#   G_sd <- round(sd(data_indeg$n), 2)
#   
#   return(c("mean" = G_avg, "median" = G_med, "sd" = G_sd))
# }
# gen.fun(g_w)


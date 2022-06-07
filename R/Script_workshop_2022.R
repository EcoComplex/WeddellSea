## Script for Workshop "Weddell Sea Dronning Maud Land"
# 14-16 June 2022

# We have decided not to focus on the comparison un-weighted / weighted properties.
# Instead, we will focus on the distribution of interaction strengths by trophic level
# and the meaning of this in the stability of the food web.


## Load packages ----

source("network_fun.r")
require(igraph)
require(multiweb)
require(tidyverse)


## Load food web database ----

# Read and check data available on GATEWAy about Weddell Sea
# Read the original database
ga <- readr::read_csv("../Data/283_2_FoodWebDataBase_2018_12_10.csv", col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
names(ga)
wedd_df <-  ga %>%filter(grepl('Weddell', foodweb.name))

# Check all species are in res.taxonomy
names(wedd_df)

# Rename and Select fields
wedd_df <- wedd_df %>% 
  rename(con_mass_mean=con.mass.mean.g., res_mass_mean=res.mass.mean.g.,
         con_taxonomy=con.taxonomy, res_taxonomy=res.taxonomy, interaction_dim=interaction.dimensionality) %>% 
  select(con_taxonomy, res_taxonomy, con_mass_mean,res_mass_mean, interaction_dim)

# Read the updated database
# Following Pawar et al. (2012) we completed missing interaction dimensionality
wedd_df <- read_delim("../Data/Wedd_mass_complete.dat", delim=" ")


## Calculate interaction intensity ----

wedd_int <- interaction_intensity(wedd_df, res_mass_mean, con_mass_mean, interaction_dim)

# Explore distribution of interaction intensity (qRC)
ggplot(wedd_int, aes(qRC)) + 
  geom_histogram(bins=50, color="darkblue", fill="white") + 
  theme_bw() + 
  scale_y_log10() 

# Convert to an igraph with weights
g <- graph_from_data_frame(wedd_int %>% 
                             select(res_taxonomy, con_taxonomy,qRC) %>% 
                             rename(weight=qRC), directed=TRUE)
E(g)$weight


# Explore interactions strength by TL ----

# Add in, out & total strengths as node (spp) attr
V(g_w)$instr <- strength(g_w, mode="in")
V(g_w)$outstr <- strength(g_w, mode="out")
V(g_w)$totalstr <- strength(g_w, mode="all")
vertex.attributes(g_w)

# Add TL as node (spp) attr
adj_mat <- as_adjacency_matrix(g_w, sparse = TRUE, attr = "weight")
tl <- round(TrophInd(as.matrix(adj_mat)), digits = 3)
V(g_w)$TL <- tl$TL
V(g_w)$Omn <- tl$OI
vertex.attributes(g_w)

# Create data frame with spp attributes
data_id <- as.data.frame(1:prop_t$Size)
data_name <- as.data.frame(V(g_w)$name)
data_totalstr <- as.data.frame(V(g_w)$totalstr)
data_instr <- as.data.frame(V(g_w)$instr)
data_outstr <- as.data.frame(V(g_w)$outstr)
data_tl <- as.data.frame(V(g_w)$TL)
data_omn <- as.data.frame(V(g_w)$Omn)
data_total <- bind_cols(data_id, data_name, data_instr, data_outstr, 
                        data_totalstr, data_tl, data_omn)
colnames(data_total) <- c("ID", "TrophicSpecies", "InStrength", 
                          "OutStrength", "TotalStrength", "TL", "Omn")

### Plot Int Strength by TL ----

# Cumulative sum & proportion of total strength
data_total <- data_total %>% 
  mutate(cumsum_str = order_by(-TotalStrength, cumsum(TotalStrength)),
         prop_str = cumsum_str/(sum(TotalStrength))) %>% 
  mutate(rank_spp = dense_rank(desc(TotalStrength)),
         prop_spp = rank_spp/nrow(data_total))

#### Total ----
(plot_totalstr_tl <- data_total %>% 
    mutate(Color = ifelse(cumsum_str < 0.6, "red", "black")) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -TL), y = TotalStrength, color = Color)) +
    geom_point() +
    scale_color_identity() +
    labs(x = "Trophic Species (descending TL)", y = "Total Strength") +
    ylim(c(0,0.018)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))

data_instr <- data_total %>% 
  mutate(cumsum_str = order_by(-InStrength, cumsum(InStrength)),
         prop_str = cumsum_str/(sum(InStrength))) %>% 
  mutate(rank_spp = dense_rank(desc(InStrength)),
         prop_spp = rank_spp/nrow(data_total))

#### In strength ----
(plot_instr_tl <- ggplot(data_total, aes(x = reorder(TrophicSpecies, -TL), y = InStrength)) +
    geom_point() +
    # geom_vline(xintercept = c(spp_nt_12+1, spp_nt_12+spp_nt_13+1, spp_nt_14+1), linetype = "longdash", colour = "red") +
    labs(x = "Trophic Species (descending TL)", y = "In Strength") +
    ylim(c(0,0.018)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))

#### Out strength ----
(plot_outstr_tl <- ggplot(data_total, aes(x = reorder(TrophicSpecies, -TL), y = OutStrength)) +
    geom_point() +
    # geom_hline(yintercept = c(spp_nt_12+1, spp_nt_12+spp_nt_13+1, spp_nt_14+1), linetype = "longdash", colour = "red") +
    labs(x = "Trophic Species (descending TL)", y = "Out Strength") +
    ylim(c(0,0.018)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))






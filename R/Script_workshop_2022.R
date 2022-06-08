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
require(NetIndices)


## Load food web database ----

# Read and check data available on GATEWAy about Weddell Sea
# Read the original database
ga <- readr::read_csv("../Data/283_2_FoodWebDataBase_2018_12_10.csv", col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
names(ga)
wedd_df <-  ga %>%
  filter(grepl('Weddell', foodweb.name))

# Rename and Select fields
wedd_df <- wedd_df %>% 
  rename(con_mass_mean = con.mass.mean.g., res_mass_mean = res.mass.mean.g.,
         con_taxonomy = con.taxonomy, res_taxonomy = res.taxonomy, interaction_dim = interaction.dimensionality) %>% 
  dplyr::select(con_taxonomy, res_taxonomy, con_mass_mean,res_mass_mean, interaction_dim)

# Read the updated database
# We completed missing interaction dimensionality following Pawar et al. (2012)
wedd_df <- read_delim("../Data/Wedd_mass_complete.dat", delim = " ")


## Calculate interaction intensity ----

wedd_int <- interaction_intensity(wedd_df, res_mass_mean, con_mass_mean, interaction_dim)

# Explore distribution of interaction intensity (qRC)
ggplot(wedd_int, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_bw() + 
  scale_y_log10() 

# Convert to an igraph with weights
g <- graph_from_data_frame(wedd_int %>% 
                             dplyr::select(res_taxonomy, con_taxonomy,qRC) %>% 
                             rename(weight=qRC), directed = TRUE)
E(g)$weight


# Explore interactions strength by TL ----

# Add in, out & total strengths as species attr
V(g)$instr <- strength(g, mode = "in")
V(g)$outstr <- strength(g, mode = "out")
V(g)$totalstr <- strength(g, mode = "all")
vertex.attributes(g)

# Add TL as species attr
adj_mat <- as_adjacency_matrix(g, sparse = TRUE, attr = "weight")
tl <- round(TrophInd(as.matrix(adj_mat)), digits = 3)
V(g)$TL <- tl$TL
V(g)$Omn <- tl$OI
vertex.attributes(g)

# Betweenness & Degree
btw <- as.data.frame(igraph::betweenness(g, directed = TRUE, weights = E(g)$weight))
deg <- as.data.frame(igraph::degree(g, mode = "total"))

# Create data frame with spp attributes
spp_name <- as.data.frame(V(g)$name)
spp_totalstr <- as.data.frame(V(g)$totalstr)
spp_instr <- as.data.frame(V(g)$instr)
spp_outstr <- as.data.frame(V(g)$outstr)
spp_tl <- as.data.frame(V(g)$TL)
spp_omn <- as.data.frame(V(g)$Omn)
spp_attr <- bind_cols(spp_name, spp_instr, spp_outstr, 
                      spp_totalstr, spp_tl, spp_omn, deg, btw)
colnames(spp_attr) <- c("TrophicSpecies", "InStrength", "OutStrength", 
                        "TotalStrength", "TL", "Omn", "Degree", "Betweenness")

### Int strength by TL ----

#### Total ----

spp_attr <- spp_attr %>% 
  mutate(cumsum_str = order_by(-TotalStrength, cumsum(TotalStrength)),
         prop_str = cumsum_str/(sum(TotalStrength))) %>% 
  mutate(rank_spp = dense_rank(desc(TotalStrength)),
         prop_spp = rank_spp/nrow(spp_attr))

(plot_totalstr_tl <- spp_attr %>% 
    mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>%
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

#### In strength ----

data_instr <- spp_attr %>% 
  mutate(cumsum_str = order_by(-InStrength, cumsum(InStrength)),
         prop_str = cumsum_str/(sum(InStrength))) %>% 
  mutate(rank_spp = dense_rank(desc(InStrength)),
         prop_spp = rank_spp/nrow(spp_attr))

(plot_instr_tl <- data_instr %>% 
    mutate(Color = ifelse(prop_str < 0.6, "red", "black")) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -TL), y = InStrength, color = Color)) +
    geom_point() +
    scale_color_identity() + 
    labs(x = "Trophic Species (descending TL)", y = "In Strength") +
    ylim(c(0,0.018)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))

#### Out strength ----

data_outstr <- spp_attr %>% 
  mutate(cumsum_str = order_by(-OutStrength, cumsum(OutStrength)),
         prop_str = cumsum_str/(sum(OutStrength))) %>% 
  mutate(rank_spp = dense_rank(desc(OutStrength)),
         prop_spp = rank_spp/nrow(spp_attr))

(plot_outstr_tl <- data_outstr %>% 
    mutate(Color = ifelse(prop_str < 0.6, "red", "black")) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -TL), y = OutStrength, color = Color)) +
    geom_point() +
    scale_color_identity() +
    labs(x = "Trophic Species (descending TL)", y = "Out Strength") +
    ylim(c(0,0.018)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))


### Total strength by Degree ----

(plot_totalstr_dg <- spp_attr %>% 
    mutate(Degree = degree(g)) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -Degree), y = TotalStrength)) +
    geom_point() +
    labs(x = "Trophic Species (descending degree)", y = "Total Strength") +
    ylim(c(0,0.018)) +
    scale_y_log10()+
    theme_bw() +
    theme(panel.grid = element_blank(),
           axis.title = element_text(size = 18, face = "bold"),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 15)))


## Plot food web ----

# By TL and degree
plot_troph_level(g, vertexSizeFactor = degree(g)*0.05)

# By TL and total interaction strength
layout_trophic <- matrix(nrow = length(V(g)), ncol = 2)
layout_trophic[, 1] <- runif(length(V(g)))
layout_trophic[, 2] <- V(g)$TL

plot.igraph(g,
            vertex.size = degree(g)*0.05,
            vertex.label = NA,
            vertex.color = V(g)$totalstr,
            layout = layout_trophic,
            edge.width = .75, edge.curved = 0.3,
            edge.arrow.size = 0.15)



# Hacer análisis

# Regresión por cuartiles

# Ajustar a power law Total Strength vs TL

# Redundancia funcional:
# Análisis multivariado con propiedades: TL, Degree, Interaction strength,
# omnivory, movement type
# Calcular betweenness pesado por Total Strength
# Diet overlap promedio





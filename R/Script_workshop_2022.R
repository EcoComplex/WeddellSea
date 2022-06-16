## Script for Workshop "Weddell Sea Dronning Maud Land"
# 14-16 June 2022

# We have decided not to focus on the comparison un-weighted / weighted properties.
# Instead, we will focus on the distribution of interaction strengths by trophic level
# and the meaning of this in the stability of the food web.


## Load packages ----

source("R/network_fun.r")
require(igraph)
require(multiweb)
require(tidyverse)
require(NetIndices)
require(plotly)

## Load food web database ----

# Read and check data available on GATEWAy about Weddell Sea
# Read the original database
ga <- readr::read_csv("Data/283_2_FoodWebDataBase_2018_12_10.csv", col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
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
wedd_df <- read_delim("Data/Wedd_mass_complete.dat", delim = " ",col_types="dccddc") %>% 
  mutate(res_mass_mean = res_mass_mean*1000,con_mass_mean = con_mass_mean*1000)


## Calculate interaction intensity ----

wedd_int <- interaction_intensity(wedd_df, res_mass_mean, con_mass_mean, interaction_dim)

# Explore distribution of interaction intensity (qRC)
ggplot(wedd_int, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  scale_y_log10()

# Convert to an igraph with weights
g <- graph_from_data_frame(wedd_int %>% 
                             dplyr::select(res_taxonomy, con_taxonomy,qRC) %>% 
                             rename(weight=qRC), directed = TRUE)
E(g)$weight


## Plot food web ----

plot_troph_level(g, vertexSizeFactor = degree(g)*0.03)


# Explore interactions strength by TL ----

# Add in, out & total strengths as species attr
V(g)$instr <- strength(g, mode = "in")
V(g)$outstr <- strength(g, mode = "out")
V(g)$totalstr <- strength(g, mode = "all")
vertex.attributes(g)

# Add weighted TL as species attr
#
adj_mat <- as_adjacency_matrix(g, sparse = TRUE, attr = "weight")

tl <- round(TrophInd(as.matrix(adj_mat)), digits = 3)
V(g)$TLw <- tl$TL
V(g)$Omnw <- tl$OI
vertex.attributes(g)


# Add unweighted TL as species attr
#
adj_mat <- as_adjacency_matrix(g, sparse = TRUE)

tl <- round(TrophInd(as.matrix(adj_mat)), digits = 3)
V(g)$TLu <- tl$TL
V(g)$Omnu <- tl$OI
vertex.attributes(g)

# Betweenness & Degree
btw <- as.data.frame(igraph::betweenness(g, directed = TRUE, weights = E(g)$weight))
deg <- as.data.frame(igraph::degree(g, mode = "total"))

# Create data frame with spp attributes
spp_name <- as.data.frame(V(g)$name)
spp_totalstr <- as.data.frame(V(g)$totalstr)
spp_instr <- as.data.frame(V(g)$instr)
spp_outstr <- as.data.frame(V(g)$outstr)
spp_tlw <- as.data.frame(V(g)$TLw)
spp_tlu <- as.data.frame(V(g)$TLu)
spp_omnw <- as.data.frame(V(g)$Omnw)
spp_attr <- bind_cols(spp_name, spp_instr, spp_outstr, 
                      spp_totalstr, spp_tlw, spp_tlu, spp_omnw, deg, btw)
colnames(spp_attr) <- c("TrophicSpecies", "InStrength", "OutStrength", 
                        "TotalStrength", "TLw", "TLu", "Omn", "Degree", "Betweenness")

### Strength by TL ----

#### Total ----

spp_attr <- spp_attr %>% 
  mutate(cumsum_str = order_by(-TotalStrength, cumsum(TotalStrength)),
         prop_str = cumsum_str/(sum(TotalStrength))) %>% 
  mutate(rank_spp = dense_rank(desc(TotalStrength)),
         prop_spp = rank_spp/nrow(spp_attr))

# Weighted trophic level
(plot_totalstr_tl <- spp_attr %>% 
    # mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -TLw), y = TotalStrength)) +  # , color = Color
    geom_point() +
    scale_y_log10() +
    scale_color_identity() +
    labs(x = "Trophic species (descending TL)", y = "Interaction strength (log scale)") +
    # ylim(c(0,0.07)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))

# Unweighted trophic level
(plot_totalstr_tl <- spp_attr %>% 
    mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -TLu), y = TotalStrength, color = Color)) +
    geom_point() +
    scale_y_log10() +
    scale_color_identity() +
    labs(x = "Trophic Species (descending TL)", y = "Interaction strength (log scale)") +
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


### Strength by Degree ----

(plot_totalstr_dg <- spp_attr %>% 
    mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>%
    ggplot(aes(x = reorder(TrophicSpecies, -Degree), y = TotalStrength, color = Color)) +
    geom_point() +
    scale_y_log10() +
    scale_color_identity() +
    labs(x = "Trophic species (descending degree)", y = "Interaction strength (log scale)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
           axis.title = element_text(size = 18, face = "bold"),
           axis.text.x = element_blank(),
           axis.text.y = element_text(size = 15)))


### Interactive plot ----

plot_3d <- plot_ly(data=spp_attr, x=~Degree, y=~TLu, 
                   z=~TotalStrength, type="scatter3d", mode="markers", 
                   color="black", marker = list(size = 5), hoverinfo = "text",
                   text = ~paste(TrophicSpecies, '<br>Degree:', Degree, '<br>Trophic level:', TLu, 
                                 '<br>Interaction strength:', round(TotalStrength,4))) %>%  
  layout(scene = list(xaxis=list(title = "Degree"),
                      yaxis=list(title = "Trophic level"),
                      zaxis=list(title = "Interaction strength", color = "red")))
plot_3d


## GAM prediction of interaction strength with TL and degree ----
#
require(mgcv)
require(gratia)

# Use weighted TL
fitw <- gam(TotalStrength ~ te(TLw, Degree, k=c(8,8)), data = spp_attr, family=tw)
plot(fitw, rug=F, pers=T, theta=45, main="Strength")
draw(fitw, residuals=T) 
appraise(fitw) 
gam.check(fitw)
summary(fitw)  # percentage of deviance explained

# Use unweighted (topological) TL
fitu <- gam(TotalStrength ~ te(TLu, Degree, k=c(8,8)), data = spp_attr, family=tw)
plot(fitu, rug=F, scheme=T, theta=45, main="Strength")
draw(fitu, residuals=T) 
appraise(fitu) 
gam.check(fitu)
summary(fitu)

# Compare GAM models (using weighted & unweighted TL)
AIC(fitw,fitu)

# Plot GAM model and observed data
# Taken from https://stackoverflow.com/questions/55047365/r-plot-gam-3d-surface-to-show-also-actual-response-values
# 2D plot

fitu_plot <- gam(TotalStrength ~ te(TLu, Degree, k=c(6,6)), data = spp_attr)  # not work with 'family=tw'
summary(fitu_plot)  # deviance % explained by the model
# Now expand it to a grid so that persp will work
steps <- 30
Degree <- with(spp_attr, seq(min(Degree), max(Degree), length = steps) )
TLu <-  with(spp_attr, seq(min(TLu), max(TLu), length = steps) )
newdat <- expand.grid(Degree = Degree, TLu = TLu)
TotalStrength <- matrix(predict(fitu_plot, newdat), steps, steps)
# Now plot it
p <- persp(Degree, TLu, TotalStrength, theta = 45, col = alpha("yellow", 0.2), 
           xlim = range(Degree), ylim = range(TLu), zlim = range(TotalStrength),
           xlab = "Degree", ylab = "Trophic level", zlab = "Interaction strength")
# To add the points, you need the same 3d transformation
obs <- with(spp_attr, trans3d(Degree, TLu, TotalStrength, p))
pred <- with(spp_attr, trans3d(Degree, TLu, fitted(fitu_plot), p))
dif <- obs[["y"]] - pred[["y"]]  # observed - predicted strength
points(obs, col = ifelse(dif < 0, alpha("blue",0.4), alpha("red",1)), pch = 16)
# Add segments to show where the points are in 3d
segments(obs$x, obs$y, pred$x, pred$y)


## To do list ----

# Quantile regression in Total Strength vs TL

# Fit power law to Total Strength vs TL

# Functional redundancy: secondary graphs

# Multivariate analysis: TL, Degree, Interaction strength, omnivory


## Figures workshop ----

# Plot food web
plot_troph_level(g, vertexSizeFactor = degree(g)*0.03)

# Plot distribution of interaction strength
ggplot(wedd_int, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  scale_y_log10()

# Calculate topological properties
calc_topological_indices(g)








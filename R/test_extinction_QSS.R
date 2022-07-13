
## Load packages
library(igraph)
library(multiweb)
library(dplyr)
library(ggplot2)
source("R/extinction_fun.R")

## Load data
load("Results/network_&_spp_attr.rda")


## Extinction sequences ----

# By trophic level
order_tl <- arrange(spp_attr_all, desc(TLu)) %>% 
  dplyr::select(TrophicSpecies, TLu)
tl_seq <- order_tl$TrophicSpecies
tl_seq <- tl_seq[1:400]

# By degree
order_deg <- arrange(spp_attr_all, desc(Degree)) %>% 
  dplyr::select(TrophicSpecies, Degree)
deg_seq <- order_deg$TrophicSpecies
deg_seq <- deg_seq[1:400]

# By mean interaction strength
order_is <- arrange(spp_attr_all, desc(AllStrength_mean)) %>% 
  dplyr::select(TrophicSpecies, AllStrength_mean)
is_seq <- order_is$TrophicSpecies
is_seq <- is_seq[1:400]


## Extinction simulations ----

# Extinctions by trophic level
tic("Tl Extinctions")
QSS_extinction_tl <- extinctions_QSS(g, tl_seq, nsim = 1000, ncores = 48, istrength = TRUE)
toc()
QSS_extinction_tl



# Extinctions by interaction strength
QSS_extinction_is <- extinctions_QSS(g, is_seq, nsim = 100, ncores = 4, istrength = TRUE)
QSS_extinction_is

# Extinctions by degree
QSS_extinction_deg <- extinctions_QSS(g, deg_seq, nsim = 100, ncores = 4, istrength = TRUE)
QSS_extinction_deg

# Extinctions topological
QSS_extinction_topol <- extinctions_QSS(g, deg_seq, nsim = 100, ncores = 4, istrength = FALSE)
QSS_extinction_topol


## Exploratory plots ----

# Results deleting by mean IS
(plot_del_QSS <- QSS_extinction_is %>% 
  mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
  ggplot(aes(x = Ext_prop, y = QSS_median)) +
  geom_line() +
  labs(x = "Proportion of deleted spp", y = "QSS median", 
       title = "Extinction sequence: decreasing mean IS"))

(plot_comp_del <- QSS_extinction_is %>% 
  mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
  ggplot(aes(x = Ext_prop, y = Components)) +
  geom_line() +
  labs(x = "Proportion of deleted spp", y = "Components", 
       title = "Extinction sequence: decreasing mean IS"))

# Results deleting by Trophic Level
(plot_del_QSS_tl <- QSS_extinction_tl %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "QSS median", 
         title = "Extinction sequence: decreasing Trophic Level"))

(plot_comp_del_tl <- QSS_extinction_tl %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = Components)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "Components", 
         title = "Extinction sequence: decreasing Trophic Level"))

# Results deleting by Trophic Level
(plot_del_QSS_deg <- QSS_extinction_deg %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "QSS median", 
         title = "Extinction sequence: decreasing Degree"))

(plot_comp_del_deg <- QSS_extinction_deg %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = Components)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "Components", 
         title = "Extinction sequence: decreasing Degree"))


# Topological results deleting by Degree
(plot_del_QSS_top <- QSS_extinction_topol %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "QSS median", 
         title = "Topological extinctions: decreasing Degree"))

(plot_comp_del_top <- QSS_extinction_topol %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = Components)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "Components", 
         title = "Topological extinctions: decreasing Degree"))


## Save data ----

save(QSS_extinction_is, QSS_extinction_tl, QSS_extinction_deg, QSS_extinction_topol,
     file = "Results/extinction_res.rda")




# Test with a small network
#
#
g1 <-   graph_from_literal( 2 -+ 1 +-3,4 -+ 1, 4-+4, 3+-3, 5-+5, 4-+6-+2, 2+-5-+3, simplify = FALSE)
c <- cluster_infomap(as.undirected(g1))
plot_troph_level(g1,vertexLabel = TRUE,vertexSizeFactor = 20,vertexSizeMin = 12, modules = TRUE, community_obj = c)
E(g1)$weight <- sample(c(.1,.2,.8,.9),gsize(g1),replace=TRUE)
tl_seq <- V(g1)$name
QSS_extinction_tl <- extinctions_QSS(g1, tl_seq, nsim = 10, ncores = 4, istrength = TRUE)

# if(count_components(g)>1){
#   dg <- components(g)
#   for(comp in unique(dg$membership)) {
#     g1 <- induced_subgraph(g, which(dg$membership == comp))
  


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
QSS_extinction_tl <- extinctions_QSS(g, tl_seq, nsim = 100, ncores = 4, istrength = TRUE)
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


## Load packages
library(igraph)
library(multiweb)
library(dplyr)
library(ggplot2)
library(tictoc)
# source("R/extinction_fun.R")

## Load data
load("Results/net_&_spp_prop_sim.rda")
load("Results/QSS_extinction_dif.rda")


## Extinction sequences ----

# By trophic level
order_tl <- arrange(spp_all_prop, desc(TL)) %>% 
  dplyr::select(TrophicSpecies, TL)
tl_seq <- order_tl$TrophicSpecies
tl_seq <- tl_seq[1:100]

# By degree
order_deg <- arrange(spp_all_prop, desc(TotalDegree)) %>% 
  dplyr::select(TrophicSpecies, TotalDegree)
deg_seq <- order_deg$TrophicSpecies
deg_seq <- deg_seq[1:100]

# By mean interaction strength
order_is <- arrange(spp_all_prop, desc(IS_mean)) %>% 
  dplyr::select(TrophicSpecies, IS_mean)
is_seq <- order_is$TrophicSpecies
is_seq <- is_seq[1:100]


## Extinction simulations ----

# Extinctions by trophic level
nsim <- 100
tic(paste("Tl Extinctions", nsim))
QSS_extinction_tl <- multiweb::calc_QSS_extinctions_seq(g, tl_seq, nsim = nsim, ncores = 48, istrength = TRUE)
toc()
QSS_extinction_tl


# Extinctions by interaction strength
QSS_extinction_is <- multiweb::calc_QSS_extinctions_seq(g, is_seq, nsim = nsim, ncores = 48, istrength = TRUE)
QSS_extinction_is

# Extinctions by degree
QSS_extinction_deg <- multiweb::calc_QSS_extinctions_seq(g, deg_seq, nsim = nsim, ncores = 48, istrength = TRUE)
QSS_extinction_deg

# Extinctions topological
# QSS_extinction_topol <- extinctions_QSS(g, deg_seq, nsim = 100, ncores = 4, istrength = FALSE)
# QSS_extinction_topol

## Save data ----

save(QSS_null_comp,QSS_null_comp_raw,QSS_extinction_dif,
     QSS_extinction_is, QSS_extinction_tl, QSS_extinction_deg, 
     file = "Results/QSS_extinction_dif.rda")

save(order_deg, order_tl, order_is, QSS_extinction_deg, QSS_extinction_tl,
     file = "Results/QSS_extinction_seq_sim.rda")

## Exploratory plots ----

# Results deleting by mean IS
#
(plot_del_QSS <- QSS_extinction_is %>% 
  mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
  ggplot(aes(x = Ext_prop, y = QSS_median)) +
  geom_line() + 
  scale_y_log10() +
  labs(x = "Proportion of deleted spp", y = "QSS median", 
       title = "Extinction sequence: by decreasing mean IS") +
   theme_classic())

(plot_comp_del <- QSS_extinction_is %>% 
  mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
  ggplot(aes(x = Ext_prop, y = Components)) +
  geom_line() +
  labs(x = "Proportion of deleted spp", y = "Components", 
       title = "Extinction sequence: by decreasing mean IS") +
  theme_classic())

# Results deleting by Trophic Level
#
(plot_del_QSS_tl <- QSS_extinction_tl %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Connectance, y = QSS_median)) +
    geom_line() +
    #labs(x = "Proportion of deleted spp", y = "QSS median", 
         #title = "Extinction sequence: by decreasing Trophic Level") +
    theme_classic())

(plot_comp_del_tl <- QSS_extinction_tl %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = Components)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "Components", 
         title = "Extinction sequence: by decreasing Trophic Level") +
    theme_classic())

# Topological results deleting by Degree
#
(plot_del_QSS_deg <- QSS_extinction_deg %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Connectance, y = QSS_median)) +
    geom_line() +
    #labs(x = "Proportion of deleted spp", y = "QSS median", 
         #title = "Extinction sequence: decreasing Degree") +
    theme_classic())

(plot_comp_del_deg <- QSS_extinction_deg %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = Components)) +
    geom_line() +
    labs(x = "Proportion of deleted spp", y = "Components", 
         title = "Extinction sequence: decreasing Degree"))







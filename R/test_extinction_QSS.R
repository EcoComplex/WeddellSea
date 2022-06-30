
## Load packages
library(igraph)
library(multiweb)
library(dplyr)
source("R/extinction_fun.R")

## Load data
load("Data/network_&_spp_attr.rda")


## Extinction sequences

# By trophic level
order_tl <- arrange(spp_attr_all, desc(TLu)) %>% 
  dplyr::select(TrophicSpecies, TLu)
tl_seq <- order_tl$TrophicSpecies
tl_seq <- tl_seq[1:10]

# By degree
order_deg <- arrange(spp_attr_all, desc(Degree)) %>% 
  dplyr::select(TrophicSpecies, Degree)
deg_seq <- order_deg$TrophicSpecies
deg_seq <- deg_seq[1:10]

# By interaction strength
order_is <- arrange(spp_attr_all, desc(AllStrength_mean)) %>% 
  dplyr::select(TrophicSpecies, AllStrength_mean)
is_seq <- order_is$TrophicSpecies
is_seq <- is_seq[1:20]


QSS_extinction_tl <- extinctions_QSS(g, tl_seq, nsim = 100, ncores = 4, istrength = TRUE)
QSS_extinction_tl

QSS_extinction_is <- extinctions_QSS(g, is_seq, nsim = 100, ncores = 4, istrength = TRUE)
QSS_extinction_is



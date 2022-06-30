

## Load packages
library(igraph)
library(multiweb)
library(dplyr)
library(ggplot2)
source("R/extinction_fun.R")

## Load data
load("Data/network_&_spp_attr.rda")


## Extinction sequences

# Sequence by trophic level
order_tl <- arrange(spp_attr_all, desc(TLu)) %>% 
  dplyr::select(TrophicSpecies, TLu)
tl_seq <- order_tl$TrophicSpecies
tl_seq <- tl_seq[1:3]

# Sequence by degree
order_deg <- arrange(spp_attr_all, desc(Degree)) %>% 
  dplyr::select(TrophicSpecies, Degree)
deg_seq <- order_deg$TrophicSpecies
deg_seq <- deg_seq[1:3]

# Sequence by interaction strength
order_is <- arrange(spp_attr_all, desc(AllStrength_mean)) %>% 
  dplyr::select(TrophicSpecies, AllStrength_mean)
deg_is <- order_is$TrophicSpecies
deg_is <- deg_is[1:3]

extinctions_QSS(g,tl_seq)

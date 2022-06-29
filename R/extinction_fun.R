#
## Function to do extinctions by Degree, Trophic Level and mean Interaction Strength
## Calculate Quasi-sign stability 'QSS' for each extinction sequence
#


## Load packages
library(igraph)
library(multiweb)
library(dplyr)
library(ggplot2)
source("R/network_fun.r")

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


# del_seq <- function(x) {

g_del <- vector(mode = "list", length = 3)
for (i in tl_seq){
  g_del <- delete_vertices(g, i)
  # Size <- vcount(g_del)
  # QSS <- multiweb::calc_QSS(g_del, nsim = 2, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
  #   mutate(median = median(QSS$maxre), Network = "All spp")
  # 
  # Data <- data.frame(Deleted = i, Size, QSS$median)
  
  }
  
  #return(Data)

#}
del_seq(g)


g_del <- vector(mode = "list", length = 3)
for(i in seq_along(g_del)) {
  for (z in seq_along(tl_seq)){
  g_del[[i]] <- delete_vertices(g, tl_seq[[z]])
  }
}

g_del

redes <- lapply(g_del, function (x) delete_vertices(g, tl_seq))
redes

g_size <- c()
for (i in 1:3){
  g_size <- c(g_size, vcount(g_del[[i]]))
}
g_size

g_QSS <- c()
for (i in 1:3){
  g_QSS <- c(g_QSS, multiweb::calc_QSS(g_del[[i]], nsim = 1, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
               mutate(median = median(g_QSS$maxre)))
}
g_QSS




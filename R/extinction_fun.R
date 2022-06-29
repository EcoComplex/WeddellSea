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


## Function

order_tl <- arrange(spp_attr_all, desc(TLu)) %>% 
  dplyr::select(TrophicSpecies, TLu)
tl_seq <- order_tl$TrophicSpecies
tl_seq <- tl_seq[1:3]


del_seq <- function(x) {
  
  for (i in tl_seq){
  
  g_del <- delete_vertices(x, i)
  Size <- vcount(g_del)
  QSS <- multiweb::calc_QSS(g_del, nsim = 2, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
    mutate(median = median(QSS$maxre), Network = "All spp")
  
  Data <- data.frame(Deleted = i, Size, QSS$median)
  
  }
  
  return(Data)

}

del_seq(g)




g_del <- delete_vertices(g, c(tl_seq))







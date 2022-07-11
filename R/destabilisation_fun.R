#
## Function to estimate destabilisation level of a species
## Calculate Quasi-sign stability 'QSS' and compare before and after extinction
## using Anderson-Darling test
#


## Load packages
library(igraph)
library(multiweb)
library(dplyr)
library(tictoc)

## Load data
load("Results/network_&_spp_attr.rda")


sp_list <- arrange(spp_attr_all, desc(AllStrength_mean)) %>% 
  dplyr::select(TrophicSpecies, AllStrength_mean)
sp_list <- sp_list$TrophicSpecies
#sp_list <- sp_list[1:2]

print("Complete Network - Nsim = 5000")
tic("QSS dif")
QSS_extinction_dif <- calc_QSS_extinction_dif(g, sp_list,ncores=48, nsim=5000, istrength = TRUE)
toc()

saveRDS(QSS_extinction_dif, "Results/QSS_extinction_dif.rds")


QSS_extinction_dif <- readRDS("Results/QSS_extinction_dif.rds")

QSS_extinction_dif %>% filter(grepl("Euphausia",Deleted))


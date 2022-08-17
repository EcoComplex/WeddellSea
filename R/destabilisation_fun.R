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
load("Results/net_&_spp_prop.rda")
load("Results/QSS_extinction_dif.rda")



sp_list <- arrange(spp_all_prop, desc(IS_mean)) %>% 
  dplyr::select(TrophicSpecies, IS_mean)
sp_list <- sp_list$TrophicSpecies

#grep("Euphausia",sp_list, value=TRUE)

#sp_list <- sp_list[c(1,87)]

nsim <- 1000
print(paste("QSS 1 sp extinction  - Nsim = ",nsim))

tic("QSS dif")
QSS_extinction_dif <- calc_QSS_extinction_dif(g, sp_list,ncores=48, nsim=nsim, istrength = TRUE)
toc()

QSS_extinction_dif <- as_tibble(QSS_extinction_dif)

save(QSS_null_comp_raw,QSS_extinction_dif,
     file = "Results/QSS_extinction_dif.rda")


QSS_extinction_dif %>% filter(grepl("Euphausia",Deleted))


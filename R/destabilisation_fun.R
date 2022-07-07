#
## Function to estimate destabilisation level of a species
## Calculate Quasi-sign stability 'QSS' and compare before and after extinction
## using Anderson-Darling test
#


## Load packages
library(igraph)
library(multiweb)
library(dplyr)

## Load data
load("Results/network_&_spp_attr.rda")


sp_list <- arrange(spp_attr_all, desc(AllStrength_mean)) %>% 
  dplyr::select(TrophicSpecies, AllStrength_mean)
sp_list <- sp_list$TrophicSpecies
sp_list <- sp_list[1:2]

calc_QSS_extinction_dif <- function(g, sp_list,nsim=10){

  comp_webs <- lapply(sp_list, function(i){
      # delete one sp and create igraph object
      g_ext <- delete_vertices(g, i)
      size <- vcount(g_ext)
      # subset mean strength for deleted sp
      # QSS for complete and deleted network
      QSS_all <- multiweb::calc_QSS(g, nsim = 3, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
        mutate(Network = "All spp")
      QSS_ext <- multiweb::calc_QSS(g_ext, nsim = 3, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
                    mutate(Network = "One sp less")
      QSS <- bind_rows(QSS_all, QSS_ext)
      # extract p-value for Anderson-Darling test comparing complete and deleted network
      ad_test <- kSamples::ad.test(maxre ~ Network, data = QSS)$ad[1,3]
      # data frame
      data.frame(Deleted = i, Network_size = size, AdTest = ad_test, QSS_all=median(QSS_all$maxre), QSS_ext=median(QSS_ext$maxre), difQSS = median(QSS_all$maxre)-median(QSS_ext$maxre))
      
  })
  bind_rows(comp_webs)
}

test <- calc_QSS_extinction_dif(g, sp_list)
sp_str <- spp_attr_all %>% 
  subset(., TrophicSpecies = i, select = AllStrength_mean)

test


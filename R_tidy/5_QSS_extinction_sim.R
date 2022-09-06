#
## Calculate Quasi-sign stability 'QSS' and compare before and after extinction simulations
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first run 'R/1_calc_interaction_strength.R' & 'R/3_species_w&uw_prop.R'


# Load packages ----

packages <- c("igraph", "multiweb", "dplyr", "tictoc")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

load("Results/net_&_spp_prop.rda")


# Run extinction simulations & estimate QSS ----

# List of species to delete
sp_list <- arrange(spp_all_prop, desc(IS_mean)) %>% 
  dplyr::select(TrophicSpecies, IS_mean)
sp_list <- sp_list$TrophicSpecies

# Be aware that this simulation might take days
# To reduce the time you can configure the number of simulations 'nsim'
nsim <- 1000
print(paste("QSS 1 sp extinction  - Nsim = ",nsim))
tic("QSS dif")
QSS_extinction_dif <- calc_QSS_extinction_dif(g, sp_list, ncores=48, nsim=nsim, istrength = TRUE)
toc()

QSS_extinction_dif <- as_tibble(QSS_extinction_dif)


# Save results ----

save(QSS_extinction_dif,
     file = "Results/QSS_extinction_dif.rda")

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

#
# Compare weighted QSS against randomized QSS with the same max weight
#
max_w <- max(E(g)$weight)
mean_w <- mean(E(g)$weight)
g1 <- g
E(g1)$weight <- max_w
mean(E(g1)$weight)

#
# qss_null
# QSS        MEing
# 0.9278 -0.002284682
if( FALSE ) {
  print("QSS Comparison mean value - Nsim = 5000")
  tic("QSS")
  qss_null <- calc_QSS(g1,ncores=48, nsim=5000, istrength = TRUE, returnRaw = FALSE)
  qss_tot  <- calc_QSS(g,ncores=48, nsim=5000, istrength = TRUE, returnRaw = FALSE)
  toc()
  
  
  QSS_null_comp <- bind_rows(QSS_null_comp, tibble(nsim=1000,QSS_null=qss_null,QSS_tot=qss_tot))
}

# QSS_null$QSS $MEing QSS_tot$QSS  $MEing
# <dbl>  <dbl>       <dbl>   <dbl>
#   1            0 0.0543           0 0.00494

print("QSS Comparison - Save Raw values - Nsim = 5000")
tic("QSS")
qss_null <- calc_QSS(g1,ncores=48, nsim=5000, istrength = TRUE, returnRaw = TRUE)
qss_tot  <- calc_QSS(g,ncores=48, nsim=5000, istrength = TRUE, returnRaw = TRUE)
toc()
QSS_null_comp_raw <- bind_rows(tibble(network="null", maxre=qss_null$maxre), 
                               tibble(network="empirical", maxre=qss_tot$maxre))

#saveRDS(QSS_null_comp, "Results/QSS_null_comp.rds")
save(QSS_null_comp,QSS_null_comp_raw,
     file = "Results/QSS_extinction_dif.rda")

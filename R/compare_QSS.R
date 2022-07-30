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
load("Results/QSS_extinction_dif.rda")

#
# Maxnull was the QSS calculated using the maximum interaction strength for the null model distribution
#
# using 2*mean interaction strength for the null gives a mean value identical to the mean interaction strength
#
QSS_null_comp_raw <- QSS_null_comp_raw %>% mutate(network=if_else(network=='null', 'maxnull', network))

#
# Compare weighted QSS against randomized QSS with the same max weight
#
max_w <- max(E(g)$weight)
mean_w <- mean(E(g)$weight)
g1 <- g
E(g1)$weight <- mean_w
g2 <- g
E(g2)$weight <- max_w

stopifnot(mean(E(g2)$weight) == max_w)

stopifnot(mean(E(g1)$weight) == mean_w)

# Testing if different number of simulations give different results
#
# qss_null
# QSS        MEing
# 0.9278 -0.002284682
if( FALSE ) {
  nsim <- 1000
  print(paste("QSS Comparison - mean values - Nsim =",nsim))
  tic("QSS")
  qss_null <- calc_QSS(g1,ncores=48, nsim=nsim, istrength = TRUE, returnRaw = FALSE)
  qss_tot  <- calc_QSS(g,ncores=48, nsim=nsim, istrength = TRUE, returnRaw = FALSE)
  toc()

  QSS_null_comp <- bind_rows(QSS_null_comp, tibble(nsim=nsim,QSS_null=qss_null,QSS_tot=qss_tot))
}



# QSS_null$QSS $MEing QSS_tot$QSS  $MEing
# <dbl>  <dbl>       <dbl>   <dbl>
#   1            0 0.0543           0 0.00494

nsim <- 1000
print(paste("QSS Comparison - Save Raw values - Nsim =",nsim))
tic("QSS")
qss_null <- calc_QSS(g1,ncores=48, nsim=nsim, istrength = TRUE, returnRaw = TRUE)
qss_null1 <- calc_QSS(g2,ncores=48, nsim=nsim, istrength = TRUE, returnRaw = TRUE)
qss_tot  <- calc_QSS(g,ncores=48, nsim=nsim, istrength = TRUE, returnRaw = TRUE)
qss_tot1  <- calc_QSS(g,ncores=48, nsim=nsim, istrength = TRUE, returnRaw = TRUE)

toc()
QSS_null_comp_raw <- bind_rows(tibble(network="meannull", maxre=qss_null$maxre), 
                               tibble(network="empirical1", maxre=qss_tot$maxre),
                               tibble(network="empirical", maxre=qss_tot1$maxre),
                               tibble(network="maxnull", maxre=qss_null1$maxre) )

#saveRDS(QSS_null_comp, "Results/QSS_null_comp.rds")
save(QSS_null_comp,QSS_null_comp_raw,
     file = "Results/QSS_extinction_dif.rda")


# Empirical 2 runs
#
emp <- QSS_null_comp_raw %>% filter(grepl("empirical",network))
ad_test <- kSamples::ad.test(maxre ~ network, data = emp)
ks.test(maxre ~ network, data = emp )

emp %>% ggplot(aes(maxre, color=network,fill=network)) + geom_density(alpha=0.5) + theme_bw()+ scale_fill_viridis_d() + scale_color_viridis_d() 
 
# 3 Networks
#
emp <- QSS_null_comp_raw %>% filter(network!="empirical")
ad_test <- kSamples::ad.test(maxre ~ network, data = emp)

# Mean Null vs Empirical
#
ad_test <- kSamples::ad.test(maxre ~ network, data = emp %>% filter(network!="maxnull"))

# Max Null vs Empirical
#
ad_test <- kSamples::ad.test(maxre ~ network, data = emp %>% filter(network!="meannull"))


## Plots
require(ggplot2)

p <- emp %>% filter(network != "maxnull") %>% 
  ggplot(aes(maxre, color=network, fill=network)) + 
  geom_density() +
#  scale_fill_viridis_d() + scale_color_viridis_d() +
  labs(x = "Maximum eigenvalue", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
p + guides(color = "none", fill = guide_legend("Network")) +
  scale_fill_discrete(labels=c('Empirical', 'Null'))

QSS_null_comp_raw %>% ggplot(aes(maxre, color=network,fill=network)) + geom_histogram(bins=30) + theme_bw() + scale_fill_viridis_d()

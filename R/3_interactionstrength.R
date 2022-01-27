## Estimate interaction strength
# Author: TIM
# Date: January 2022


## Load packages ----

pkg_list <- c("dplyr", "readr", "igraph", "multiweb", "NetIndices", "ggplot2")

lapply(pkg_list, FUN = function(pkg_list) {
  do.call("require", list(pkg_list)) 
})


## Load data ----

# Load total (data_total) & curated data (weddell_df) from 1_datawrangling.R
load("../Data/Weddell_database.rda")
load("../Data/g_unweighted_att.rda")


# Resource/Consumer ratio
weddell_ratio <- data_total %>% 
  mutate(mass.ratio = res.mass.mean.g./con.mass.mean.g.)


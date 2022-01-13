## Data wrangling
# Author: TIM
# Date: January 2022


## Load packages ----

pkg_list <- c("dplyr", "readr")

lapply(pkg_list, FUN = function(pkg_list) {
  do.call("require", list(pkg_list)) 
})

## Load data from GATEWAy Database ----





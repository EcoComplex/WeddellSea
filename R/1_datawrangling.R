## Data wrangling
# Author: TIM
# Date: January 2022


## Load packages ----

pkg_list <- c("dplyr", "readr", "igraph")

lapply(pkg_list, FUN = function(pkg_list) {
  do.call("require", list(pkg_list)) 
})


## Load data from GATEWAy Database ----

data_total <- readr::read_csv("../Data/FoodWebDataBase_GATEWAy.csv", col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
str(data_total)
spec(data_total)  # explore columns

# Filter Weddell Sea food web
unique(data_total$foodweb.name)
grep('Weddell', unique(data_total$foodweb.name), value=TRUE)

weddell_df <-  data_total %>%
  filter(grepl('Weddell', foodweb.name))


## Save curated data ----

save(weddell_df, data_total, file = "../Data/Weddell_database.rda")

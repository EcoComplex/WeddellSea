#
## Summary interaction strength IS by sp
#

# Interactions list: data frame
wedd_int

#
# Using 'wedd_int' data frame
#

# Consumers (incoming interactions)
con_int <- wedd_int %>% 
  dplyr::select(con_taxonomy, qRC) %>% 
  dplyr::rename(TrophicSpecies = con_taxonomy, IS = qRC)

# Resources (outgoing interactions)
res_int <- wedd_int %>% 
  dplyr::select(res_taxonomy, qRC) %>% 
  dplyr::rename(TrophicSpecies = res_taxonomy, IS = qRC)

# Cannibalistic interactions
canib_int <- wedd_int %>% 
  dplyr::filter(duplicated(cbind(con_taxonomy, res_taxonomy)))
canib_count <- wedd_int %>% 
  group_by(con_taxonomy, res_taxonomy) %>% 
  mutate(Canibal = n() > 1)

is_simple(g)
any_loop(g)
loops <- length(which_multiple(g) = TRUE)  # not run

# All interactions
all_int <- bind_rows(con_int, res_int)
total_int <- bind_rows(con_int, res_int) %>% 
  dplyr::group_by(TrophicSpecies) %>% 
  dplyr::summarize(AllStrength_mean = mean(IS),
            AllStrength_Q1 = quantile(IS, 0.25),
            AllStrength_Q3 = quantile(IS, 0.75),
            AllStrength_count = n())

# Add to spp attributes data frame 'spp_all'
total_int
spp_all

spp_attr_all <- spp_all %>% 
  left_join(total_int) %>% 
  dplyr::select(TrophicSpecies, Degree, Betweenness, TLu, Omnu, meanTrophicSimil,
                InStrength_sum, OutStrength_sum, AllStrength_sum,
                AllStrength_mean, AllStrength_Q1, AllStrength_Q3, AllStrength_count)

save(wedd_df, wedd_int, g, all_int, spp_attr_all,
     file = "Results/network_&_spp_attr.rda")








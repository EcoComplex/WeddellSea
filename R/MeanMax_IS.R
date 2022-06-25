#
## Mean and max interaction strength IS by sp
#

# Interactions list: data frame
wedd_int
# Food web: igraph object
g

m <- as_adjacency_matrix(g, attr = "weight", names = TRUE)
m <- as.matrix(m)

prey_names <- rownames(m)
pred_names <- colnames(m)

m_by_row <- as.data.frame(m[prey_names, ])
head(m_by_row)
m_by_col <- as.data.frame(m[, pred_names])
head(m_by_col)



#
# Using 'wedd_int' data frame
#

# Consumers (incoming interactions)
con <- wedd_int %>% 
  dplyr::select(con_taxonomy, qRC) %>% 
  dplyr::group_by(con_taxonomy) %>% 
  mutate(mean_str_con = mean(qRC), max_str_con = max(qRC)) %>% 
  dplyr::distinct(con_taxonomy, .keep_all = TRUE) %>% 
  dplyr::rename(Species = con_taxonomy)

# Resources (outgoing interactions)
res <- wedd_int %>% 
  dplyr::select(res_taxonomy, qRC) %>% 
  dplyr::group_by(res_taxonomy) %>%
  mutate(mean_str_res = mean(qRC), max_str_res = max(qRC)) %>%
  dplyr::distinct(res_taxonomy, .keep_all = TRUE) %>% 
  dplyr::rename(Species = res_taxonomy)

all_spp <- bind_rows(con, res) %>% 
  group_by(Species) %>% 
  mutate(mean_final = mean(mean_str_con, mean_str_res),
         max_final = max(max_str_con, max_str_res))

# all_spp[is.na(all_spp)] <- -1




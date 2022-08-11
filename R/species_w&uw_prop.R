#
## Weighted (interaction strength) & unweighted species properties
#


# Load packages ----
packages <- c("dplyr", "igraph", "NetIndices", "cheddar")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Weighted properties ----

## Load data ----
load("Results/interaction_estimation.rda")


## Add interaction strength per species ----

# Select interaction strength for Consumers (incoming)
con_int <- wedd_int_pd %>% 
  dplyr::select(con_taxonomy, qRC) %>% 
  dplyr::rename(TrophicSpecies = con_taxonomy, IS = qRC)

# Select interaction strength for Resources (outcoming)
res_int <- wedd_int_pd %>% 
  dplyr::select(res_taxonomy, qRC) %>% 
  dplyr::rename(TrophicSpecies = res_taxonomy, IS = qRC)

# Cannibalistic interactions
# canib_int <- wedd_int %>% 
#   dplyr::filter(duplicated(cbind(con_taxonomy, res_taxonomy)))
# canib_count <- wedd_int %>% 
#   group_by(con_taxonomy, res_taxonomy) %>% 
#   mutate(Canibal = n() > 1)
# 
# is_simple(g)
# any_loop(g)
# loops <- length(which_multiple(g) = TRUE)  # not run

# Bind interaction strengths for each species
# Calculate mean, Q1 & Q3
all_int <- bind_rows(con_int, res_int)
total_int <- bind_rows(con_int, res_int) %>% 
  dplyr::group_by(TrophicSpecies) %>% 
  dplyr::summarize(IS_mean = mean(IS),
                   IS_median = median(IS),
                   IS_Q1 = quantile(IS, 0.25),
                   IS_Q3 = quantile(IS, 0.75),
                   IS_sum_tot = sum(IS),
                   IS_max = max(IS),
            Check_NumbInt = n())


# Unweighted properties ----

## Load data ----

# Convert interaction list to an igraph with weights
g <- graph_from_data_frame(wedd_int_pd %>% 
                             dplyr::select(res_taxonomy, con_taxonomy, qRC) %>% 
                             rename(weight=qRC), directed = TRUE)
E(g)$weight
# Add sum of in & out interaction strengths per species
V(g)$IS_sum_in <- strength(g, mode = "in")
V(g)$IS_sum_out <- strength(g, mode = "out")
vertex_attr_names(g)


## Calculate properties ----

# Trophic level
adj_mat <- as_adjacency_matrix(g, sparse = TRUE)
tl <- round(TrophInd(as.matrix(adj_mat)), digits = 3)
V(g)$TL <- tl$TL
V(g)$Omn <- tl$OI
vertex.attributes(g)

# Degree
V(g)$Deg <- degree(g, mode = "total")
V(g)$InDeg <- degree(g, mode = "in")
V(g)$OutDeg <- degree(g, mode = "out")
vertex_attr_names(g)

# Trophic similarity
source("R/igraph_cheddar.R")  # load function to convert igraph to cheddar object
igraph_to_cheddar(g)

cc <- LoadCommunity("Community")
ts <- TrophicSimilarity(cc)
mts <- tibble(TrophicSpecies = rownames(ts), meanTrophicSimil = colMeans(ts))  # data frame


# All spp properties ----

# Weighted
spp_IS_sum_in <- data.frame(V(g)$name, V(g)$IS_sum_in)
colnames(spp_IS_sum_in) <- c("TrophicSpecies", "IS_sum_in")
spp_IS_sum_out <- data.frame(V(g)$name, V(g)$IS_sum_out)
colnames(spp_IS_sum_out) <- c("TrophicSpecies", "IS_sum_out")
spp_w_prop <- total_int %>% 
  left_join(spp_IS_sum_in) %>% 
  left_join(spp_IS_sum_out)

# Unweighted
spp_name <- as.data.frame(V(g)$name)
spp_totdeg <- as.data.frame(V(g)$Deg)
spp_indeg <- as.data.frame(V(g)$InDeg)
spp_outdeg <- as.data.frame(V(g)$OutDeg)
spp_tl <- as.data.frame(V(g)$TL)
spp_omn <- as.data.frame(V(g)$Omn)
spp_tropsimil <- mts[,2]  

spp_uw_prop <- bind_cols(spp_name, spp_totdeg, spp_indeg, spp_outdeg, spp_tl,
                         spp_omn, spp_tropsimil)
colnames(spp_uw_prop) <- c("TrophicSpecies", "TotalDegree", "InDegree", "OutDegree",
                        "TL", "Omn", "meanTrophicSimil")

# Bind weighted & unweighted properties
spp_w_prop
spp_uw_prop

spp_all_prop <- spp_uw_prop %>% 
  left_join(spp_w_prop) %>% 
  dplyr::select(TrophicSpecies, TotalDegree, InDegree, OutDegree, TL, Omn, meanTrophicSimil,
                IS_mean, IS_median, IS_max, IS_Q1, IS_Q3, IS_sum_tot, IS_sum_in, IS_sum_out)


# Save results ----

save(g, spp_w_prop, spp_uw_prop,
     file = "Results/net_&_spp_prop.rda")

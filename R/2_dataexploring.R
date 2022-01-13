## Data exploring
# Author: TIM
# Date: January 2022


## Load packages ----

pkg_list <- c("dplyr", "readr", "igraph", "multiweb", "NetIndices")

lapply(pkg_list, FUN = function(pkg_list) {
  do.call("require", list(pkg_list)) 
})


## Load data ----

# Load total (data_total) & curated data (weddell_df) from 1_datawrangling.R
load("../Data/Weddell_database.rda")


## Create igraph object for unweighted network ----

edgelist_unw <- as.matrix(weddell_df %>%
                     dplyr::select(res.taxonomy, con.taxonomy))
g_unw  <- graph_from_edgelist(edgelist_unw, directed  = T)
nodelist_unw <- as.data.frame(vertex_attr(g_unw))  # node list


## Calculate topological (unw) indices ----

# Common network indices
multiweb::calc_topological_indices(g_unw)

# Degree
indeg <- degree(g_unw, mode = "in")
outdeg <- degree(g_unw, mode = "out")
totdeg <- degree(g_unw, mode = "total")
V(g_unw)$indegree <- indeg
V(g_unw)$outdegree <- outdeg
V(g_unw)$totdegree <- totdeg

# Species type (basal, intermediate, top)
basal <- V(g_unw)[indegree == 0]
V(g_unw)[basal]$sp.type <- "basal"
top <- V(g_unw)[outdegree == 0]
V(g_unw)[top]$sp.type <- "top"
int <- V(g_unw)[indegree != 0 & outdegree != 0]
V(g_unw)[int]$sp.type <- "int"

# Trophic level & omnivory
tl_omn <-TrophInd(get.adjacency(g_unw, sparse = F))
V(g_unw)$trophic.level <- tl_omn$TL
V(g_unw)$omn.ind <- tl_omn$OI

vertex.attributes(g_unw)

# Generality
gen.fun <- function(g){
  pred <- degree(g, mode = "in") > 0
  G <- round(sum(degree(g, mode = "in")[pred] / sum(pred)), 2)
  return(G)
}
gen.fun(g_unw)

# Vulnerability
vul.fun <- function(g){
  prey <- degree(g, mode = "out") > 0
  V <- round(sum(degree(g, mode = "out")[prey]) / sum(prey), 2)
  return(V)
}
vul.fun(g_unw)


## Save unweighted igraph w/attributes

save(g_unw, file = "../Data/g_unweighted_att.rda")

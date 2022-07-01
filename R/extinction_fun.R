#
## Function to do extinctions by Degree, Trophic Level and mean Interaction Strength
## Calculate Quasi-sign stability 'QSS' for each extinction sequence
#


#' QSS for an extinction sequence
#' 
#' @param g_del igraph object for the deletion sequence
#' @param seq   vector or list with the nodes for the extinction sequence, the 
#'              extinctions are incremental. 
#' @param nsim        number of simulations used in QSS calculation
#' @param ncores      number of cores used to calculate QSS
#' @param istrength   parameter istrength for QSS function 
#'
#' @return data frame with the median QSS, the number of nodes and the name
#'         of the last deleted node
#' @export
#'
#' @examples

extinctions_QSS <- function(g_del, seq, nsim=100, ncores=0, istrength=FALSE){
  
  qdel <- lapply(seq, function(i){
    g_del <<- delete_vertices(g_del, i) 
    size <- vcount(g_del)
    con <- multiweb::calc_topological_indices(g_del)$Connectance
    comp <- multiweb::calc_topological_indices(g_del)$Components
    QSS <- multiweb::calc_QSS(g_del, nsim = nsim, ncores = ncores, istrength = istrength, returnRaw = TRUE) %>% summarize(QSS_median = median(maxre), QSS_prop = sum(maxre <0)/nsim) %>%  
      mutate(Size = size, Connectance = con, Components = comp, Last_deleted = i)
  })
  qdel <- bind_rows(qdel)
}



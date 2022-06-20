#
# Transform g object ('igraph') into community files ('cheddar')
#

require(igraph)
require(cheddar)
require(readr)

igraph_to_cheddar <- function(g) {
  
  # create a folder "Community" in working directory to store needed files
  dir.create("Community")
  # obtain trophic links file
  el <- as_edgelist(g)
  readr::write_delim(data.frame(resource = el[,1], consumer = el[,2]), file = "Community/trophic.links.csv", delim = ",")
  # obtain nodes file
  n <- vertex_attr(g, "name")
  readr::write_delim(data.frame(node = n), file = "Community/nodes.csv", delim = ",")
  # create 'food web name' file
  p <- data.frame(title = "Food Web name")
  readr::write_delim(p, file = "Community/properties.csv", delim = ",")
  # load and create Community object
  cc <<- LoadCommunity("Community")

}

igraph_to_cheddar(g)

# Calculate trophic similarity
ts <- TrophicSimilarity(cc)
# Generate a data frame with mean trophic similarity for each species
mts <- tibble(TrophicSpecies=rownames(ts), meanTrophicSimil=colMeans(ts))

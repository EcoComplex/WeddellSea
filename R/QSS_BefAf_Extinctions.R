#
## Quasi-sign stability 'QSS' before & after extinctions
#


## Load packages
library(igraph)
library(multiweb)


## Load data
load("Data/network_&_spp_attr.rda")


## Weddell Sea FW before extinction
g

# Create FWs without key species according to literature
g_eup_sup <- delete_vertices(g, "Euphausia superba")
g_orc_orc <- delete_vertices(g, "Orcinus orca")
g_ple_ant <- delete_vertices(g, "Pleuragramma antarcticum")
g_mac_hol <- delete_vertices(g, "Macrourus holotrachys")

# QSS before
QSS <- multiweb::calc_QSS(g, nsim = 100, ncores = 4, istrength = TRUE) %>% 
  mutate(Network = "All spp") %>% rename(QSS=QSS, MEing=MEing)
# QSS without Euphausia superba
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_eup_sup, nsim = 100, ncores = 4, istrength = TRUE) %>% 
              mutate(Network = "Without Eup. sup.") %>% rename(QSS = QSS, MEing=MEing))
# QSS without Orcinus orca
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_orc_orc, nsim = 100, ncores = 4, istrength = TRUE) %>% 
              mutate(Network = "Without Orc. orc.") %>% rename(QSS = QSS, MEing=MEing))
# QSS without Pleuragramma antarcticum
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_ple_ant, nsim = 100, ncores = 4, istrength = TRUE) %>% 
              mutate(Network = "Without Ple. ant.") %>% rename(QSS = QSS, MEing=MEing))
# QSS without Macrourus holotrachys
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_mac_hol, nsim = 100, ncores = 4, istrength = TRUE) %>% 
              mutate(Network = "Without Mac. hol.") %>% rename(QSS = QSS, MEing=MEing))

# Compare properties with & without extinctions
websTbl_all <- calcTopologicalIndices(g) %>% 
  mutate(Network = "All spp") %>% 
  dplyr::select(Network, everything())
websTbl_ext1 <- calcTopologicalIndices(g_eup_sup) %>% 
  mutate(Network = "Without Eup. sup.") %>% 
  dplyr::select(Network, everything())
websTbl_ext2 <- calcTopologicalIndices(g_orc_orc) %>% 
  mutate(Network = "Without Orc. orc.") %>% 
  dplyr::select(Network, everything())
websTbl_ext3 <- calcTopologicalIndices(g_ple_ant) %>% 
  mutate(Network = "Without Ple. ant.") %>% 
  dplyr::select(Network, everything())
websTbl_ext4 <- calcTopologicalIndices(g_mac_hol) %>% 
  mutate(Network = "Without Mac. hol.") %>% 
  dplyr::select(Network, everything())

websTbl <- bind_rows(websTbl_all, websTbl_ext1, websTbl_ext2, 
                     websTbl_ext3, websTbl_ext4)

websTbl <- websTbl %>% 
  inner_join(QSS) %>% 
  dplyr::select(Network, MEing, QSS, everything())













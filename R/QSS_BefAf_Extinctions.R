#
## Quasi-sign stability 'QSS' before & after extinctions
#


## Load packages
library(igraph)
library(multiweb)
library(dplyr)
library(ggplot2)

## Load data
load("Data/network_&_spp_attr.rda")


## Weighted ----
## Weddell Sea FW before extinctions
g

# Create FWs without key species according to literature
g_eup_sup <- delete_vertices(g, "Euphausia superba")
g_orc_orc <- delete_vertices(g, "Orcinus orca")
g_ple_ant <- delete_vertices(g, "Pleuragramma antarcticum")
g_mac_hol <- delete_vertices(g, "Macrourus holotrachys")
g_phyto <- delete_vertices(g, "Phytodetritus")
g_mes_ham <- delete_vertices(g, "Mesonychoteuthis hamiltoni")

# QSS all spp
QSS <- multiweb::calc_QSS(g, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
  mutate(Network = "All spp")
# QSS without Euphausia superba
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_eup_sup, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Eup. sup."))
# QSS without Orcinus orca
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_orc_orc, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Orc. orc."))
# QSS without Pleuragramma antarcticum
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_ple_ant, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Ple. ant."))
# QSS without Macrourus holotrachys
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_mac_hol, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Mac. hol."))
# QSS without Macrourus holotrachys
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_phyto, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Phytodet."))
# QSS without Mesonychoteuthis hamiltoni
QSS <- QSS %>% 
  bind_rows(multiweb::calc_QSS(g_mes_ham, nsim = 1000, ncores = 4, istrength = TRUE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Mes. ham."))

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
                     websTbl_ext3)

websTbl <- websTbl %>% 
  inner_join(QSS) %>% 
  dplyr::select(Network, MEing, QSS, everything())

# Plots
ggplot(QSS, aes(x = maxre, color = Network)) +
  geom_density()
ggplot(QSS, aes(x = maxre, color = Network)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ Network)

# Statistical analysis
ks.test(maxre ~ Network, data = QSS)
kSamples::ad.test(maxre ~ Network, data = QSS)


## Un-weighted ----

# QSS all spp
QSS_u <- multiweb::calc_QSS(g, nsim = 1000, ncores = 4, istrength = FALSE, returnRaw = TRUE) %>% 
  mutate(Network = "All spp")
# QSS without Euphausia superba
QSS_u <- QSS_u %>% 
  bind_rows(multiweb::calc_QSS(g_eup_sup, nsim = 1000, ncores = 4, istrength = FALSE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Eup. sup."))
# QSS without Orcinus orca
QSS_u <- QSS_u %>% 
  bind_rows(multiweb::calc_QSS(g_orc_orc, nsim = 1000, ncores = 4, istrength = FALSE, returnRaw = TRUE) %>% 
              mutate(Network = "Without Orc. orc."))

# Plots
ggplot(QSS_u, aes(x = maxre, color = Network)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ Network)

ggplot(QSS_u, aes(x = maxre, color = Network)) +
  geom_density()

# Statistical analysis
kSamples::ad.test(maxre ~ Network, data = QSS_u)







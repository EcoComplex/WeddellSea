#
## Extinction sequences
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022-2023
#


# Load packages -----------------------------------------------------------
packages <- c("igraph", "multiweb", "dplyr", "ggplot2", "tictoc", "tidyverse")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ---------------------------------------------------------------
load("Results/net_&_spp_prop_sim.rda")
load("Results/QSS_extinction_seq_sim.rda")
#load("Results/QSS_extinction_dif.rda")


# Extinction simulations ----------------------------------------------------

## Trophic level ----
# Sequence
order_tl <- arrange(spp_all_prop, desc(TL)) %>% 
  dplyr::select(TrophicSpecies, TL)
tl_seq <- order_tl$TrophicSpecies
tl_seq <- tl_seq[1:490]
# Simulation
nsim <- 100
tic(paste("Tl Extinctions", nsim))
QSS_extinction_tl <- multiweb::calc_QSS_extinctions_seq(g, tl_seq, nsim = nsim, ncores = 48, istrength = TRUE)
toc()

## Degree ----
# Sequence
order_deg <- arrange(spp_all_prop, desc(TotalDegree)) %>% 
  dplyr::select(TrophicSpecies, TotalDegree)
deg_seq <- order_deg$TrophicSpecies
deg_seq <- deg_seq[1:490]
# Simulation
QSS_extinction_deg <- multiweb::calc_QSS_extinctions_seq(g, deg_seq, nsim = nsim, ncores = 48, istrength = TRUE)

## Interaction strength ----
# Sequence
order_is <- arrange(spp_all_prop, desc(IS_mean)) %>% 
  dplyr::select(TrophicSpecies, IS_mean)
is_seq <- order_is$TrophicSpecies
is_seq <- is_seq[1:490]
# Simulation
QSS_extinction_is <- multiweb::calc_QSS_extinctions_seq(g, is_seq, nsim = nsim, ncores = 48, istrength = TRUE)


# Exploratory plots -------------------------------------------------------
# Proportion of deleted spp vs QSS
(plot_Prop_QSS <- QSS_extinction_is %>% 
   mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
   ggplot(aes(x = Ext_prop, y = QSS_median)) +
   geom_line() + 
   scale_y_log10() +
   labs(x = "Proportion of deleted spp", y = "QSS median") +
   theme_classic())
# Connectance vs QSS
(plot_Conn_QSS <- QSS_extinction_deg %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Connectance, y = QSS_median)) +
    geom_line() +
    theme_classic())

(plot_Conn_QSS <- QSS_extinction_deg %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() +
    theme_classic())

(plot_Conn_QSS <- QSS_extinction_is %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() +
    theme_classic())

(plot_Conn_QSS <- QSS_extinction_tl %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() +
    theme_classic())
plotly::ggplotly(plot_Conn_QSS)


## Functional Group ----

## Read full taxonomic Classification  
classall_df <- readRDS("Data/WeddellSea_clasification.rds")

### Porifera & Ascidians ----
por_asc <- classall_df %>% filter(phylum == "Porifera" | class =="Ascidiacea") %>% dplyr::select(Taxon) %>% deframe()
# Simulation
nsim <- 100
QSS_extinction_grp <- multiweb::calc_QSS_extinctions_seq(g, por_asc, nsim=nsim, ncores=8, istrength = TRUE)
# Plot
(plot_Conn_QSS <- QSS_extinction_grp %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() + geom_smooth(method=lm) +
    labs(title="Porifera & Ascidians deleted (n=46)", x="Proportion of deleted spp") +
    theme_classic())
plotly::ggplotly(plot_Conn_QSS)
# Difference with original network
nsim <- 1000
por_asc_dif <- data.frame()
por_asc_dif <- bind_rows(por_asc_dif,
                     calc_QSS_extinction_dif_grp(g, por_asc, nsim, ncores=8, istrength=TRUE) %>% mutate(deleted_grp="Pori - Asci"))

### Mammals & Sea birds ----
av_mm <- classall_df %>% filter(class == "Aves" | class =="Mammalia") %>% dplyr::select(Taxon) %>% deframe()
# Simulation
nsim <- 100
QSS_extinction_grp_am <- multiweb::calc_QSS_extinctions_seq(g, av_mm, nsim=nsim, ncores=8, istrength = TRUE)
# Plot
(plot_Conn_QSS <- QSS_extinction_grp_am %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() + geom_smooth(method=lm) +
    labs(title="Mammals & Sea birds deleted (n=30)", x="Proportion of deleted spp") +
    theme_classic())
plotly::ggplotly(plot_Conn_QSS)
# Difference with original network
nsim <- 1000
av_mm_dif <- data.frame()
av_mm_dif <- bind_rows(av_mm_dif,
                     calc_QSS_extinction_dif_grp(g, av_mm, nsim, ncores=8, istrength=TRUE) %>% mutate(deleted_grp="Mamm - Birds"))

### Myctopids & Euphasids ----
m_e <- classall_df %>% filter(family == "Myctophidae" | family =="Euphausiidae") %>% dplyr::select(Taxon) %>% deframe()
# Simulation
nsim <- 100
QSS_extinction_grp_me <- multiweb::calc_QSS_extinctions_seq(g, m_e, nsim=nsim, ncores=8, istrength = TRUE)
# Plot
(plot_Conn_QSS <- QSS_extinction_grp_me %>% 
    mutate(Network_prop = Size/490, Ext_prop = (490-Size)/490) %>% 
    ggplot(aes(x = Ext_prop, y = QSS_median)) +
    geom_line() + geom_smooth(method=lm) +
    labs(title="Myctopids & Euphasids deleted (n=10)", x="Proportion of deleted spp") +
    theme_classic())
plotly::ggplotly(plot_Conn_QSS)
# Difference with original network
nsim <- 1000
m_e_dif <- data.frame()
m_e_dif <- bind_rows(m_e_dif,
                     calc_QSS_extinction_dif_grp(g, m_e, nsim, ncores=8, istrength=TRUE) %>% mutate(deleted_grp="Mycto - Eupha"))

# All data
dif_all <- bind_rows(por_asc_dif, av_mm_dif, m_e_dif)

# Save results ------------------------------------------------------------

# save(QSS_null_comp,QSS_null_comp_raw,QSS_extinction_dif,
#      QSS_extinction_is, QSS_extinction_tl, QSS_extinction_deg, 
#      file = "Results/QSS_extinction_dif.rda")

save(order_deg, order_tl, order_is, QSS_extinction_deg, QSS_extinction_tl, QSS_extinction_is,
     QSS_extinction_grp, QSS_extinction_grp_am, QSS_extinction_grp_me, QSS_extinction_grp_dif,
     dif_all, file = "Results/QSS_extinction_seq_sim.rda")


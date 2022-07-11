
library(ggplot2)
library(dplyr)

## Load data
data <- readRDS(file = "Results/QSS_extinction_dif.rds")
load("Results/network_&_spp_attr.rda")


## Habitat data
hab_data <- read.csv(file = "Data/WeddellSeaHabitat.csv")
hab_data_inc <- hab_data %>% 
  dplyr::select(species, Environment) %>% 
  rename(TrophicSpecies = species, Habitat = Environment)
hab_data_comp <- hab_data_inc %>% 
  add_row(TrophicSpecies = "Phytodetritus", Habitat = "Benthic") %>%        # missing species
  add_row(TrophicSpecies = "Sediment", Habitat = "Benthic") %>%             # missing species
  add_row(TrophicSpecies = "Iphimediella  cyclogena", Habitat = "Benthic")  # missing species


## Join extinction results (QSS), Habitat & spp properties
all_data <- data %>% 
  rename(TrophicSpecies = Deleted) %>% 
  left_join(spp_attr_all) %>% 
  left_join(hab_data_comp)


## Correlation p-values (Kolmogorov vs Anderson Darling)
ggplot(data, aes(x = KS_pvalue, y = Ad_pvalue)) +
  geom_point()


## Relationship btw QSS difference and spp properties
ggplot(all_data, aes(x = TLu, y = difQSS)) +
  geom_point()

ggplot(all_data, aes(x = Degree, y = difQSS)) +
  geom_point() +
  scale_x_log10()

ggplot(all_data, aes(x = AllStrength_mean, y = difQSS)) +
  geom_point() +
  scale_x_log10()

ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS)) +
  geom_point()

ggplot(all_data, aes(x = Omnu, y = difQSS)) +
  geom_point() +
  scale_x_log10()

ggplot(all_data, aes(x = Habitat, y = difQSS)) +
  geom_boxplot()


## 20 + 20 influencers

# Destabilisers by presence
des_spp_pres <- all_data %>% 
  arrange(desc(difQSS)) %>% 
  slice(1:20) %>% 
  mutate(Destabilization = "Presence")
# Destabilisers by abscence
des_spp_abs <- all_data %>% 
  arrange(desc(-difQSS)) %>% 
  slice(1:20) %>% 
  mutate(Destabilization = "Abscence")

inf_spp <- bind_rows(des_spp_pres, des_spp_abs)


ggplot(inf_spp, aes(x = TLu, y = difQSS)) +
  geom_point()

ggplot(inf_spp, aes(x = Habitat, y = difQSS)) +
  geom_point()

ggplot(inf_spp, aes(x = meanTrophicSimil, y = difQSS)) +
  geom_point()

ggplot(inf_spp, aes(x = Degree, y = difQSS)) +
  geom_point()

ggplot(inf_spp, aes(x = AllStrength_mean, y = difQSS)) +
  geom_point() +
  scale_x_log10()


ggplot(wedd_int, aes(qRC)) +
  geom_histogram()

inf_spp_name <- inf_spp$TrophicSpecies
inf_spp_dist <- all_int[all_int$TrophicSpecies %in% inf_spp_name, ] %>% 
  left_join(inf_spp)
unique(inf_spp_dist$TrophicSpecies)  

library(ggjoy)
p <- ggplot(inf_spp_dist, aes(x = IS, y = reorder(TrophicSpecies, TLu, decreasing=T))) +
  geom_density_ridges(aes(fill = Destabilization), jittered_points = TRUE) +
  scale_x_log10()
p + labs(fill = "Destabilization")

ggplot(inf_spp, aes(y = TLu)) +
  geom_boxplot() +
  geom_jitter()


## GAM

library(mgcv)
library(gratia)

all_fit <- gam(difQSS ~ te(TLu, Degree, k=c(8,8)), data = all_data)
plot(all_fit, rug=F, pers=T, theta=45, main="Strength")
draw(all_fit, residuals=T) 
appraise(all_fit) 
gam.check(all_fit)
summary(all_fit)  # percentage of deviance explained






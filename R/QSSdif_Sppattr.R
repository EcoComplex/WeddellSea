#
## Relationship btw QSS difference and spp properties
# Positive values indicate spp that are destabilizers
# Negative values indicate spp that are stabilizers


#
## Load packages

library(dplyr)
library(ggplot2)


#
## Load data

load("Results/QSS_extinction_dif.rda")
load("Results/network_&_spp_attr.rda")


#
## Compare QSS empirical and null

p <- ggplot(QSS_null_comp_raw, aes(x = maxre, fill = network, color = network)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  labs(x = "QSS eigenvalues", y = "Frequency", fill = "Network") +
  theme_classic()
p + guides(color = FALSE)


#
## Merge QSS results with spp attributes

# Habitat data
hab_data <- read.csv(file = "Data/WeddellSeaHabitat.csv")
hab_data_inc <- hab_data %>% 
  dplyr::select(species, Environment) %>% 
  rename(TrophicSpecies = species, Habitat = Environment)
hab_data_comp <- hab_data_inc %>% 
  add_row(TrophicSpecies = "Phytodetritus", Habitat = "Benthic") %>%            # missing species
  add_row(TrophicSpecies = "Sediment", Habitat = "Benthic") %>%                 # missing species
  add_row(TrophicSpecies = "Iphimediella  cyclogena", Habitat = "Benthic") %>%  # missing species
  mutate(Habitat = case_when(Habitat == "Bathydemersal" ~ "Demersal", TRUE ~ Habitat))  # replace 'Bathydemersal' with 'Demersal'

# Join extinction results (QSS), Habitat & spp properties
QSS_data <- QSS_extinction_dif

all_data <- QSS_data %>% 
  rename(TrophicSpecies = Deleted) %>% 
  left_join(spp_attr_all) %>% 
  left_join(hab_data_comp)


#
## By trophic level

ggplot(all_data, aes(x = TLu, y = difQSS)) +
  geom_point(aes(color = ifelse(cluster == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  labs(color = "Group", shape = "Stability impact", x = "Trophic level", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

#
## By degree

ggplot(all_data, aes(x = Degree, y = difQSS)) +
  geom_point(aes(color = ifelse(cluster == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_x_log10() +
  labs(color = "Group", shape = "Stability impact", x = "Degree (log scale)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

#
## By mean interaction strength

ggplot(all_data, aes(x = log(AllStrength_mean), y = difQSS)) +
  geom_point(aes(color = ifelse(cluster == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  labs(color = "Group", shape = "Stability impact", x = "log(mean Interaction strength)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

#
## By trophic similarity

ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS)) +
  geom_point(aes(color = ifelse(cluster == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  labs(color = "Group", shape = "Stability impact", x = "Trophic similarity", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

#
## By omnivory

ggplot(all_data, aes(x = Omnu, y = difQSS)) +
  geom_point(aes(color = ifelse(cluster == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  labs(color = "Group", shape = "Stability impact", x = "Omnivory", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

#
## By habitat

ggplot(all_data, aes(x = Habitat, y = difQSS)) +
  geom_violin(fill = "gray") +
  geom_point(aes(color = ifelse(cluster == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  #scale_color_manual(values = c("black", "blue")) +
  labs(color = "Group", shape = "Stability impact", x = "Habitat", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15))


#
## Statistical correlations 

# p-values (Kolmogorov vs Anderson Darling)
ggplot(all_data, aes(x = KS_pvalue, y = Ad_pvalue)) +
  geom_point() +
  labs(x = "Kolmogorov-Smirnov p-value", y = "Anderson-Darling p-value") +
  theme_classic()

# Subset sp with QSS difference significant (AD test p-value < 0.01)
key_sp <- all_data %>% 
  filter(., Ad_pvalue < 0.01) %>% 
  dplyr::select(TrophicSpecies, AllStrength_mean, TLu, Degree, meanTrophicSimil,
                Habitat, difQSS, Ad_pvalue) %>% 
  arrange(Ad_pvalue)

# p-values and QSS difference
library(reshape2)
comp_QSS_pvalue <- melt(all_data, measure.vars = c("Ad_pvalue", "KS_pvalue"))
ggplot(comp_QSS_pvalue, aes(x = value, y = abs(difQSS), color = variable)) +
  geom_point() +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("#E69F00","#56B4E9"), labels = c('Anderson','Kolmogorov')) +
  labs(x = "p-value", y = "abs(QSS difference)", color = "Test") +
  theme_classic()


#
## Relationship btw QSS difference and spp properties
# Positive values indicate spp that are destabilizers
# Negative values indicate spp that are stabilizers


#
# Load packages ----

library(dplyr)
library(ggplot2)


#
# Load data ----

load("Results/QSS_extinction_dif.rda")
load("Results/net_&_spp_prop.rda")
load("Results/cluster_data.rds")


#
# Compare QSS empirical-null ----

p <- ggplot(QSS_null_comp_raw, aes(x = maxre, fill = network, color = network)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  labs(x = "QSS eigenvalues", y = "Frequency", fill = "Network") +
  theme_classic()
p + guides(color = "none")


#
# Merge QSS w/ spp prop ----

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

# Join extinction results (QSS), Habitat, clustering & spp properties
QSS_data <- QSS_extinction_dif

# Rename clusters by relative IS
cluster_data["cluster_mean"][cluster_data["cluster_mean"] == "1"] <- "High"
cluster_data["cluster_mean"][cluster_data["cluster_mean"] == "2"] <- "Low"

cluster_data["cluster_median"][cluster_data["cluster_median"] == "1"] <- "High"
cluster_data["cluster_median"][cluster_data["cluster_median"] == "3"] <- "Medium"
cluster_data["cluster_median"][cluster_data["cluster_median"] == "2"] <- "Low"

cluster_data["cluster_sum_tot"][cluster_data["cluster_sum_tot"] == "1"] <- "High"
cluster_data["cluster_sum_tot"][cluster_data["cluster_sum_tot"] == "2"] <- "Low"

cluster_data["cluster_max"][cluster_data["cluster_max"] == "1"] <- "Single"

all_data <- QSS_data %>% 
  rename(TrophicSpecies = Deleted) %>% 
  left_join(hab_data_comp) %>% 
  left_join(cluster_data) %>% 
  left_join(spp_all_prop)


# QSS vs spp prop ----

## By mean interaction strength ----

IS_QSS <- ggplot(all_data, aes(x = log(IS_mean), y = difQSS)) +
  geom_point(aes(color = cluster_mean,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Group", shape = "Stability impact", x = "log(mean Interaction strength)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
IS_QSS

# For IS_max
ggplot(all_data, aes(x = log(IS_max), y = difQSS)) +
  geom_point(aes(color = cluster_max,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(color = "Group", shape = "Stability impact", x = "log(max Interaction strength)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

## By trophic level ----

TL_QSS <- ggplot(all_data, aes(x = TL, y = difQSS)) +
  geom_point(aes(color = cluster_mean,   # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Group", shape = "Stability impact", x = "Trophic level", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
TL_QSS

# For IS_max
ggplot(all_data, aes(x = TL, y = difQSS)) +
  geom_point(aes(color = cluster_max,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(color = "Group", shape = "Stability impact", x = "Trophic level", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

## By degree ----

DEG_QSS <- ggplot(all_data, aes(x = TotalDegree, y = difQSS)) +
  geom_point(aes(color = cluster_mean,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  scale_x_log10() +
  labs(color = "Group", shape = "Stability impact", x = "Degree (log scale)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
DEG_QSS

# For IS_max
ggplot(all_data, aes(x = TotalDegree, y = difQSS)) +
  geom_point(aes(color = cluster_max,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  scale_x_log10() +
  labs(color = "Group", shape = "Stability impact", x = "Degree (log scale)", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

## By trophic similarity ----

TS_QSS <- ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS)) +
  geom_point(aes(color = cluster_mean, 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Group", shape = "Stability impact", x = "Trophic similarity", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
TS_QSS

# For IS_max
ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS)) +
  geom_point(aes(color = cluster_max,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(color = "Group", shape = "Stability impact", x = "Trophic similarity", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

## By omnivory ----

ggplot(all_data, aes(x = Omn, y = difQSS)) +
  geom_point(aes(color = ifelse(cluster_mean == "High", "High IS", "Low IS"), 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(color = "Group", shape = "Stability impact", x = "Omnivory", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

## By habitat ----

HAB_QSS <- ggplot(all_data, aes(x = Habitat, y = difQSS)) +
  geom_violin(fill = "grey90") +
  geom_point(aes(color = cluster_mean,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  labs(color = "Group", shape = "Stability impact", x = "Habitat", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 12))
HAB_QSS

# For IS_max
ggplot(all_data, aes(x = Habitat, y = difQSS)) +
  geom_violin(fill = "gray") +
  geom_point(aes(color = cluster_max,  # 'cluster_median' or 'cluster_sum_tot'
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(color = "Group", shape = "Stability impact", x = "Habitat", y = "Stability difference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15))

## Figure 6 ----
library(ggpubr)
library(grid)
IS_QSS
TL_QSS
DEG_QSS
TS_QSS
HAB_QSS

legend <- ggplot(all_data, aes(x = meanTrophicSimil, y = difQSS))+
  geom_point(aes(color = cluster_mean, 
                 shape = ifelse(Ad_pvalue < 0.01, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_shape_manual(values = c(19, 2), labels = c("Non-significant", "Significant")) +
  lims(x = c(0,0), y = c(0,0)) +
  theme_void() +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  10),
        legend.title = element_text(size = 12, face = "bold")) +
  guides(color = guide_legend(title = "Group"), shape = guide_legend(title = "Stability impact"))

Fig6_QSS <- ggarrange(IS_QSS + rremove("ylab"), TL_QSS + rremove("ylab"), 
                      DEG_QSS + rremove("ylab"), TS_QSS + rremove("ylab"), 
                      HAB_QSS + rremove("ylab"), legend,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol=3, nrow=2)
Fig6_QSS <- annotate_figure(Fig6_QSS, left = text_grob("Stability difference", rot = 90, 
                                          vjust = 1, size = 18))

ggsave(filename = "Manuscript/Fig.6_QSSDif.png", plot = Fig6_QSS, 
       width = 10, units = "in", dpi = 600, bg = "white")


# Statistical correlations ----

# p-values (Kolmogorov vs Anderson Darling)
ggplot(all_data, aes(x = KS_pvalue, y = Ad_pvalue)) +
  geom_point() +
  labs(x = "Kolmogorov-Smirnov p-value", y = "Anderson-Darling p-value") +
  theme_classic()

# Subset sp with QSS difference significant (AD test p-value < 0.01)
key_sp <- all_data %>% 
  filter(., Ad_pvalue < 0.01) %>% 
  left_join(cluster_data) %>% 
  dplyr::select(TrophicSpecies, IS_mean, TL, TotalDegree, meanTrophicSimil,
                Habitat, difQSS, Ad_pvalue) %>% 
  arrange(Ad_pvalue)


# p-value - QSS difference ----

library(reshape2)
comp_QSS_pvalue <- melt(all_data, measure.vars = c("Ad_pvalue", "KS_pvalue"))
ggplot(comp_QSS_pvalue, aes(x = value, y = abs(difQSS), color = variable)) +
  geom_point() +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("#E69F00","#56B4E9"), labels = c('Anderson','Kolmogorov')) +
  labs(x = "p-value", y = "abs(QSS difference)", color = "Test") +
  theme_classic()


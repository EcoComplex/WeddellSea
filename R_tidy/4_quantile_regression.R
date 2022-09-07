#
## Quantile regressions for relationships btw interaction strength and spp properties
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first run 'R/1_calc_interaction_strength.R' & 'R/3_species_w&uw_prop.R'


# Load packages ----

packages <- c("ggplot2", "dplyr", "lsmeans", "olsrr", "factoextra", "cluster")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

load("Results/net_&_spp_prop.rda")


# Distribution of species' interaction strength ----

# Distribution of log IS
ggplot(spp_all_prop, aes(x = log(IS_mean))) + geom_density() + theme_bw()
ggplot(spp_all_prop, aes(x = log(IS_mean))) + geom_histogram(bins=50) + theme_bw()


# Cluster by interaction strength ----

# Determine and visualize the optimal number of clusters using total within sum of square
data.scaled <- scale(log(spp_all_prop$IS_mean))
fviz_nbclust(data.scaled, kmeans, method = "wss")
# Calculate gap statistic based on number of clusters
gap_stat <- clusGap(data.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
# Plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
# Perform k-means clustering with k = predicted clusters (gap_stat)
km <- kmeans(data.scaled, centers = 2, nstart = 25)

# Add cluster assignment to original data
cluster_data <- cbind(spp_all_prop, cluster_mean = km$cluster)

# Plot clusters
ggplot(cluster_data, aes(y = log(IS_mean), x=cluster_mean, color=cluster_mean)) + 
  geom_jitter() + theme_bw() + scale_color_viridis_c()

p <- ggplot(cluster_data, aes(x = log(IS_mean))) + 
  geom_density(aes(group=as.factor(cluster_mean), color=as.factor(cluster_mean), fill=as.factor(cluster_mean)), alpha=0.3) + 
  labs(x = "log(mean Interaction strength)", y = "Frequency", fill = "Cluster") +
  theme_bw() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
p + guides(color = "none")

# Rename clusters by relative interaction strength
cluster_data["cluster_mean"][cluster_data["cluster_mean"] == "1"] <- "High"
cluster_data["cluster_mean"][cluster_data["cluster_mean"] == "2"] <- "Low"

# Save cluster results
save(cluster_data, file = "Results/cluster_data.rds")

# Regression by clusters (aka groups) ----

## Trophic level ----

cl_IS_TL <- ggplot(cluster_data, aes(x = TL, y = log(IS_mean), color = cluster_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Trophic level", y = "log(mean Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
cl_IS_TL

# Test regression significance by group
TL_lm <- lm(log(IS_mean) ~ TL * cluster_mean, data = cluster_data)
summary(TL_lm)
TL_eq <- lstrends(TL_lm, "cluster_mean", var="TL")
TL_eq
pairs(TL_eq)  # significance btw group slopes
ols_test_normality(TL_lm)  # check normality of residuals
ols_plot_resid_qq(TL_lm)  # Q-Q plot


## Degree ----

cl_IS_DEG <- ggplot(cluster_data, aes(x = TotalDegree, y = log(IS_mean), color = cluster_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Degree (log scale)", y = "log(mean Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
cl_IS_DEG

# Test regression significance by group
DEG_lm <- lm(log(IS_mean) ~ TotalDegree * cluster_mean, data = cluster_data)
summary(DEG_lm)
DEG_eq <- lstrends(DEG_lm, "cluster_mean", var="TotalDegree")
DEG_eq
pairs(DEG_eq)  # significance btw group slopes
ols_test_normality(DEG_lm)  # check normality of residuals
ols_plot_resid_qq(DEG_lm)  # Q-Q plot


## Trophic similarity ----

cl_IS_TS <- ggplot(cluster_data, aes(x = meanTrophicSimil, y = log(IS_mean), color = cluster_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Trophic similarity", y = "log(mean Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
cl_IS_TS

# Test regression significance by group
TS_lm <- lm(log(IS_mean) ~ meanTrophicSimil * cluster_mean, data = cluster_data)
summary(TS_lm)
TS_eq <- lstrends(TS_lm, "cluster_mean", var="meanTrophicSimil")
TS_eq
pairs(TS_eq)  # significance btw group slopes
ols_test_normality(TS_lm)  # check normality of residuals
ols_plot_resid_qq(TS_lm)  # Q-Q plot


## Habitat ----

cl_IS_HAB <- ggplot(cluster_data, aes(x = Habitat, y = log(IS_mean))) +
  geom_violin(fill = "grey70", alpha = 0.5) +
  geom_point(shape=21, aes(fill = factor(cluster_mean))) +
  scale_fill_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Habitat", y = "log(mean Interaction Strength)", fill = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12))
cl_IS_HAB


# Save plots ----

save(cl_IS_TL, cl_IS_DEG, cl_IS_TS, cl_IS_HAB,
     file = "Results/single_plots.rda")

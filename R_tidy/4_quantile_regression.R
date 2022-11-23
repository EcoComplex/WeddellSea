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


# Check cluster by interaction strength ----

# Determine and visualize the optimal number of clusters using total within sum of square
data.scaled <- scale(log(spp_all_prop$IS_mean))
fviz_nbclust(data.scaled, kmeans, method = "wss")
# Calculate gap statistic based on number of clusters
gap_stat <- clusGap(data.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
# Plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
# Perform k-means clustering with k = predicted clusters (gap_stat)
km <- kmeans(data.scaled, centers = 1, nstart = 25)
# NO CLUSTERING FOUND! #


# Linear Regression ----

## Trophic level ----

cl_IS_TL <- ggplot(spp_all_prop, aes(x = TL, y = log(IS_mean))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Trophic level", y = "log(mean Interaction Strength)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
cl_IS_TL

# Test regression significance
TL_lm <- lm(log(IS_mean) ~ TL, data = spp_all_prop)
summary(TL_lm)
ols_test_normality(TL_lm)  # check normality of residuals
ols_plot_resid_qq(TL_lm)  # Q-Q plot


## Degree ----

cl_IS_DEG <- ggplot(spp_all_prop, aes(x = TotalDegree, y = log(IS_mean))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  labs(x = "Degree (log scale)", y = "log(mean Interaction Strength)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
cl_IS_DEG

# Test regression significance
DEG_lm <- lm(log(IS_mean) ~ TotalDegree, data = spp_all_prop)
summary(DEG_lm)
ols_test_normality(DEG_lm)  # check normality of residuals
ols_plot_resid_qq(DEG_lm)  # Q-Q plot


## Trophic similarity ----

cl_IS_TS <- ggplot(spp_all_prop, aes(x = meanTrophicSimil, y = log(IS_mean))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Trophic similarity", y = "log(mean Interaction Strength)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
cl_IS_TS

# Test regression significance
TS_lm <- lm(log(IS_mean) ~ meanTrophicSimil, data = spp_all_prop)
summary(TS_lm)
ols_test_normality(TS_lm)  # check normality of residuals
ols_plot_resid_qq(TS_lm)  # Q-Q plot


## Habitat ----

cl_IS_HAB <- ggplot(spp_all_prop, aes(x = Habitat, y = log(IS_mean))) +
  geom_violin(fill = "grey70", alpha = 0.5) +
  geom_point(shape=19) +
  labs(x = "Habitat", y = "log(mean Interaction Strength)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12))
cl_IS_HAB


# Save plots ----

save(cl_IS_TL, cl_IS_DEG, cl_IS_TS, cl_IS_HAB,
     file = "Results/single_plots_sep22.rda")

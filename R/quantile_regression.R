#
## Quantile regression for relationships btw interaction strength and spp properties
#


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


# Species IS distribution ----

# Distribution of log IS
ggplot(spp_all_prop, aes(x = log(IS_mean))) + geom_density() + theme_bw() 
ggplot(spp_all_prop, aes(x = log(IS_mean))) + geom_histogram(bins=50) + theme_bw()  # mean
ggplot(spp_all_prop, aes(x = log(IS_median))) + geom_histogram(bins=50) + theme_bw()  # median
ggplot(spp_all_prop, aes(x = log(IS_max))) + geom_histogram(bins=50) + theme_bw()  # max
ggplot(spp_all_prop, aes(x = log(IS_sum_tot))) + geom_histogram(bins=50) + theme_bw()  # total sum
ggplot(spp_all_prop, aes(x = log(IS_sum_in))) + geom_histogram(bins=50) + theme_bw()  # in sum
ggplot(spp_all_prop, aes(x = log(IS_sum_out))) + geom_histogram(bins=50) + theme_bw()  # out sum


# Clustering ----

# Explore clustering by IS (K-means & Gap statistic)
# by 'IS_mean', 'IS_median', 'IS_max', 'IS_sum_tot'

data.scaled <- scale(log(spp_all_prop$IS_max))
fviz_nbclust(data.scaled, kmeans, method = "wss")
# calculate gap statistic based on number of clusters
gap_stat <- clusGap(data.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
# plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)
# perform k-means clustering with k = predicted clusters (gap_stat)
km <- kmeans(data.scaled, centers = 1, nstart = 25)

# add cluster assignment to original data
cluster_data <- cbind(spp_all_prop, cluster_mean = km$cluster)

cluster_data$cluster_median <- km$cluster
cluster_data$cluster_sum_tot <- km$cluster
cluster_data$cluster_max <- km$cluster

# Save cluster results
save(cluster_data, file = "Results/cluster_data.rds")


# Load cluster data

load("Results/cluster_data.rds")


# Plot groups
ggplot(cluster_data, aes(y = log(IS_max), x=cluster_max, color=cluster_max)) + 
  geom_jitter() + theme_bw() + scale_color_manual(values = c("#3a5e8cFF"))  # + scale_color_viridis_c()

p <- ggplot(cluster_data, aes(x = log(IS_max))) + 
  geom_density(aes(group=as.factor(cluster_max), color=as.factor(cluster_max), fill=as.factor(cluster_max)), alpha=0.3) + 
  labs(x = "log(max Interaction strength)", y = "Frequency", fill = "Cluster") +
  theme_bw() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
p + guides(color = "none")


# Explore relationships ----
# btw interaction strength (mean, median, sum_tot) & spp properties

attach(cluster_data)
datatable <- data.frame(TotalDegree, TL, Omn, meanTrophicSimil, 
                        log(IS_mean), log(IS_median), log(IS_sum_tot), log(IS_max))
cor(datatable, method = "spearman")
pairs(datatable, col="blue", main="Scatterplots")


## Trophic level ----

# Test
# QR <- rq(log(AllStrength_mean) ~ TLu, tau=c(0.15, 0.85), data=spp_attr_all)
# summary(QR, se= "boot")
# QR_15 <- rq(log(AllStrength_mean) ~ TLu, tau=0.15, data=spp_attr_all)
# QR_85 <- rq(log(AllStrength_mean) ~ TLu, tau=0.85, data=spp_attr_all)
# anova(QR_15, QR_85)  # test difference btw quantiles 15 & 85

# Plot
qr_IS_TL <- ggplot(cluster_data, aes(x = TL, y = log(IS_sum_tot))) +
  geom_point(shape=21, aes(fill = factor(cluster_sum_tot))) +
  scale_fill_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  #geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Trophic level", y = "log(sum Interaction Strength)", fill = "Group") +  # , color = "Quantiles" 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_TL


## Degree ----

# Test
# QRD <- rq(log(AllStrength_mean) ~ Degree, tau=c(0.15, 0.85), data=spp_attr_all)
# summary(QRD, se= "boot")
# QRD_15 <- rq(log(AllStrength_mean) ~ Degree, tau=0.15, data=spp_attr_all)
# QRD_85 <- rq(log(AllStrength_mean) ~ Degree, tau=0.85, data=spp_attr_all)
# anova(QRD_15, QRD_85)

# Plot
qr_IS_DEG <- ggplot(cluster_data, aes(x = TotalDegree, y = log(IS_sum_tot))) +
  geom_point(shape=21, aes(fill = factor(cluster_sum_tot))) +
  scale_fill_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  scale_x_log10() +
  #geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Degree (log scale)", y = "log(sum Interaction Strength)", fill = "Group") +  # , color = "Quantiles"
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_DEG


## Trophic Simil ----

# Test
# QRTS <- rq(log(AllStrength_mean) ~ meanTrophicSimil, tau=c(0.15,0.5, 0.85), data=spp_attr_all)
# summary(QRTS, se= "boot")
# QRS_15 <- rq(log(AllStrength_mean) ~ meanTrophicSimil, tau=0.15, data=spp_attr_all)
# QRS_85 <- rq(log(AllStrength_mean) ~ meanTrophicSimil, tau=0.85, data=spp_attr_all)
# anova(QRS_15, QRS_85)

# Plot
qr_IS_TS <- ggplot(cluster_data, aes(x = meanTrophicSimil, y = log(IS_sum_tot))) +
  geom_point(shape=21, aes(fill = factor(cluster_sum_tot))) +
  scale_fill_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  #geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Mean Trophic Similarity", y = "log(sum Interaction Strength)", fill = "Group") +  # , color = "Quantiles"
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_TS


## Habitat ----

ggplot(cluster_data, aes(x = Habitat, y = log(IS_sum_tot))) +
  geom_violin(fill = "grey70", alpha = 0.5) +
  geom_point(shape=21, aes(fill = factor(cluster_sum_tot))) +
  scale_fill_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Habitat", y = "log(sum Interaction Strength)", fill = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15))

# For IS_max
ggplot(cluster_data, aes(x = Habitat, y = log(IS_max), color = cluster_max)) +
  geom_violin(fill = "grey70", alpha = 0.5) +
  geom_point() +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(x = "Habitat", y = "log(max Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15))


# Regression by groups ----

# Rename clusters by relative IS
cluster_data["cluster_mean"][cluster_data["cluster_mean"] == "1"] <- "High"  # 'cluster_sum_tot'
cluster_data["cluster_mean"][cluster_data["cluster_mean"] == "2"] <- "Low"  # 'cluster_sum_tot'

cluster_data["cluster_median"][cluster_data["cluster_median"] == "1"] <- "High"
cluster_data["cluster_median"][cluster_data["cluster_median"] == "3"] <- "Medium"
cluster_data["cluster_median"][cluster_data["cluster_median"] == "2"] <- "Low"

cluster_data["cluster_max"][cluster_data["cluster_max"] == "1"] <- "Cluster"

## Trophic level ----

cl_IS_TL <- ggplot(cluster_data, aes(x = TL, y = log(IS_sum_tot), color = cluster_sum_tot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Trophic level", y = "log(sum Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
cl_IS_TL

# For IS_max
ggplot(cluster_data, aes(x = TL, y = log(IS_max), color = cluster_max)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(x = "Trophic level", y = "log(max Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

# Test regression significance by group
TL_lm <- lm(log(IS_sum_tot) ~ TL * cluster_sum_tot, data = cluster_data)
summary(TL_lm)
TL_eq <- lstrends(TL_lm, "cluster_sum_tot", var="TL")
TL_eq
pairs(TL_eq)  # significance btw group slopes
ols_test_normality(TL_lm)  # check normality of residuals
ols_plot_resid_qq(TL_lm)  # Q-Q plot


## Degree ----

cl_IS_DEG <- ggplot(cluster_data, aes(x = TotalDegree, y = log(IS_sum_tot), color = cluster_sum_tot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Degree (log scale)", y = "log(sum Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
cl_IS_DEG

# For IS_max
ggplot(cluster_data, aes(x = TotalDegree, y = log(IS_max), color = cluster_max)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(x = "Degree (log scale)", y = "log(max Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

# Test regression significance by group
DEG_lm <- lm(log(IS_sum_tot) ~ TotalDegree * cluster_sum_tot, data = cluster_data)
summary(DEG_lm)
DEG_eq <- lstrends(DEG_lm, "cluster_sum_tot", var="TotalDegree")
DEG_eq
pairs(DEG_eq)  # significance btw group slopes
ols_test_normality(DEG_lm)  # check normality of residuals
ols_plot_resid_qq(DEG_lm)  # Q-Q plot


## Trophic similarity ----

cl_IS_TS <- ggplot(cluster_data, aes(x = meanTrophicSimil, y = log(IS_sum_tot), color = cluster_sum_tot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("#541352FF", "#ffcf20FF"), labels = c("High IS", "Low IS")) +
  labs(x = "Trophic similarity", y = "log(sum Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
cl_IS_TS

# For IS_max
ggplot(cluster_data, aes(x = meanTrophicSimil, y = log(IS_max), color = cluster_max)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("#3a5e8cFF"), labels = c("Single")) +
  labs(x = "Trophic similarity", y = "log(max Interaction Strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

# Test regression significance by group
TS_lm <- lm(log(IS_sum_tot) ~ meanTrophicSimil * cluster_sum_tot, data = cluster_data)
summary(TS_lm)
TS_eq <- lstrends(TS_lm, "cluster_sum_tot", var="meanTrophicSimil")
TS_eq
pairs(TS_eq)  # significance btw group slopes
ols_test_normality(TS_lm)  # check normality of residuals
ols_plot_resid_qq(TS_lm)  # Q-Q plot


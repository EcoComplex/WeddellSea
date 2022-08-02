#
## Quantile regression for relationships btw interaction strength and spp attributes
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

load("Results/network_&_spp_attr.rda")


# Separate spp by Int str ----

# Distribution of log IS
#
ggplot(spp_attr_all, aes(x = log(AllStrength_mean))) + geom_density() + theme_bw() 
ggplot(spp_attr_all, aes(x = log(AllStrength_mean))) + geom_histogram(bins=50) + theme_bw() 

#
# Kmeans to separate the two groups
#
#km <- kmeans(log(spp_attr_all$AllStrength_mean), 2)
#spp_attr_all$cluster <- km$cluster

data.scaled <- scale(log(spp_attr_all$AllStrength_mean))
fviz_nbclust(data.scaled, kmeans, method = "wss")
gap_stat <- clusGap(data.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)


# Plot 2 groups
#ggplot(spp_attr_all, aes(y = log(AllStrength_mean),x=cluster,color=cluster)) + geom_jitter() + theme_bw() + scale_color_viridis_c()
ggplot(spp_attr_all, aes(x = log(AllStrength_mean), color=cluster)) + 
  geom_density() + 
  labs(x = "log(mean Interaction strength)", y = "Density", color = "Group") +
  theme_bw() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))


# Explore relationships ----
# btw interaction strength & spp attributes

attach(spp_attr_all)
datatable <- data.frame(Degree, TLu, Omnu, meanTrophicSimil, 
                        AllStrength_sum, AllStrength_mean)
cor(datatable)
pairs(datatable, col="blue", main="Scatterplots")


## Trophic level ----

# Test
QR <- rq(log(AllStrength_mean) ~ TLu, tau=c(0.15, 0.85), data=spp_attr_all)
summary(QR, se= "boot")
QR_15 <- rq(log(AllStrength_mean) ~ TLu, tau=0.15, data=spp_attr_all)
QR_85 <- rq(log(AllStrength_mean) ~ TLu, tau=0.85, data=spp_attr_all)
anova(QR_15, QR_85)  # test difference btw quantiles 15 & 85

# Plot
qr_IS_TL <- ggplot(spp_attr_all, aes(x = TLu, y = log(AllStrength_mean))) +
  geom_point(shape=21, aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("red", "blue"), labels = c("High IS", "Low IS")) +
  geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Trophic level", y = "mean Interaction Strength (log scale)", color = "Quantiles", fill = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_TL


## Degree ----

# Test
QRD <- rq(log(AllStrength_mean) ~ Degree, tau=c(0.15, 0.85), data=spp_attr_all)
summary(QRD, se= "boot")
QRD_15 <- rq(log(AllStrength_mean) ~ Degree, tau=0.15, data=spp_attr_all)
QRD_85 <- rq(log(AllStrength_mean) ~ Degree, tau=0.85, data=spp_attr_all)
anova(QRD_15, QRD_85)

# Plot
qr_IS_DEG <- ggplot(spp_attr_all, aes(x = Degree, y = log(AllStrength_mean))) +
  geom_point(shape=21, aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("red", "blue"), labels = c("High IS", "Low IS")) +
  geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Degree", y = "Interaction strength (log scale)", color = "Quantiles", fill = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_DEG


## Trophic Simil ----

# Test
QRTS <- rq(log(AllStrength_mean) ~ meanTrophicSimil, tau=c(0.15,0.5, 0.85), data=spp_attr_all)
summary(QRTS, se= "boot")
QRS_15 <- rq(log(AllStrength_mean) ~ meanTrophicSimil, tau=0.15, data=spp_attr_all)
QRS_85 <- rq(log(AllStrength_mean) ~ meanTrophicSimil, tau=0.85, data=spp_attr_all)
anova(QRS_15, QRS_85)

# Plot
qr_IS_TS <- ggplot(spp_attr_all, aes(x = meanTrophicSimil, y = log(AllStrength_mean))) +
  geom_point(shape=21, aes(fill = factor(cluster))) +
  scale_fill_manual(values = c("red", "blue"), labels = c("High IS", "Low IS")) +
  geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Mean Trophic Similarity", y = "Interaction strength (log scale)", color = "Quantiles", fill = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_TS


# Regression by groups ----

# Rename clusters
#spp_attr_all["cluster"][spp_attr_all["cluster"] == "2"] <- "High"
#spp_attr_all["cluster"][spp_attr_all["cluster"] == "1"] <- "Low"


## Trophic level ----

cl_IS_TL <- ggplot(spp_attr_all, aes(x = TLu, y = log(AllStrength_mean), color = cluster)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("High IS", "Low IS")) +
  labs(x = "Trophic level", y = "log(mean Interaction strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
cl_IS_TL

# Test regression significance by group
TL_lm <- lm(log(AllStrength_mean) ~ TLu * cluster, data = spp_attr_all)
summary(TL_lm)
TL_eq <- lstrends(TL_lm, "cluster", var="TLu")
TL_eq
pairs(TL_eq)  # significance btw group slopes
ols_test_normality(TL_lm)  # check normality of residuals
ols_plot_resid_qq(TL_lm)  # Q-Q plot


## Degree ----

cl_IS_DEG <- ggplot(spp_attr_all, aes(x = Degree, y = log(AllStrength_mean), color = cluster)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("High IS", "Low IS")) +
  labs(x = "Degree", y = "log(mean Interaction strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
cl_IS_DEG

# Test regression significance by group
DEG_lm <- lm(log(AllStrength_mean) ~ Degree * cluster, data = spp_attr_all)
summary(DEG_lm)
DEG_eq <- lstrends(DEG_lm, "cluster", var="Degree")
DEG_eq
pairs(DEG_eq)  # significance btw group slopes
ols_test_normality(DEG_lm)  # check normality of residuals
ols_plot_resid_qq(DEG_lm)  # Q-Q plot


## Trophic similarity ----

cl_IS_TS <- ggplot(spp_attr_all, aes(x = meanTrophicSimil, y = log(AllStrength_mean), color = cluster)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_discrete(labels = c("High IS", "Low IS")) +
  labs(x = "Trophic similarity", y = "log(mean Interaction strength)", color = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
cl_IS_TS

# Test regression significance by group
TS_lm <- lm(log(AllStrength_mean) ~ meanTrophicSimil * cluster, data = spp_attr_all)
summary(TS_lm)
TS_eq <- lstrends(TS_lm, "cluster", var="meanTrophicSimil")
TS_eq
pairs(TS_eq)  # significance btw group slopes
ols_test_normality(TS_lm)  # check normality of residuals
ols_plot_resid_qq(TS_lm)  # Q-Q plot



# Save results ----

save(all_int, g, spp_attr_all, wedd_df, wedd_int, 
     file="Results/network_&_spp_attr.rda")

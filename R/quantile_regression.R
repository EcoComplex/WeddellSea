#
## Quantile regression for relationships btw interaction strength and spp attributes
#


# Load packages ----

packages <- c("SparseM", "quantreg", "ggplot2")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load data ----

load("Results/network_&_spp_attr.rda")


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
  geom_point() + 
  geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Trophic level", y = "mean Interaction Strength (log scale)", color = "Quantiles") +
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
  geom_point() +
  geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Degree", y = "Interaction strength (log scale)", color = "Quantiles") +
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
  geom_point() +
  geom_quantile(quantiles = c(0.15, 0.85), size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  labs(x = "Mean Trophic Similitude", y = "Interaction strength (log scale)", color = "Quantiles") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_TS


#
# Distribution of log IS
#
ggplot(spp_attr_all, aes( x = log(AllStrength_mean))) + geom_density() + theme_bw() 
ggplot(spp_attr_all, aes( x = log(AllStrength_mean))) + geom_histogram(bins=50) + theme_bw() 

#
# Kmeans to separate the two groups
#
km <- kmeans(log(spp_attr_all$AllStrength_mean),2)

spp_attr_all$cluster <- km$cluster

ggplot(spp_attr_all, aes( y = log(AllStrength_mean),x=cluster,color=cluster)) + geom_jitter() + theme_bw() + scale_color_viridis_c()
ggplot(spp_attr_all, aes( x = log(AllStrength_mean),color=factor(cluster))) + geom_density() + theme_bw() + scale_color_viridis_d()

save(all_int,g,spp_attr_all,wedd_df,wedd_int, file="Results/network_&_spp_attr.rda")

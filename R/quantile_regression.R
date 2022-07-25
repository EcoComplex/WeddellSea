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

IntStr <- cbind(AllStrength_mean)
TL <- cbind(TLu)

hist(IntStr, prob=TRUE, col = "blue", border = "black", xlab = "mean Interaction Strength")
lines(density(IntStr))

# Lineal regression
OLSreg <- lm(IntStr ~ TL)
summary(OLSreg)

# Quantile 25 & 75
QR_2575 <- rq(IntStr ~ TL, tau=c(0.25, 0.75))
sumQR_2575_boot <- summary(QR_2575, se = "boot")  # test significance of regression

QR_25 <- rq(IntStr ~ TL, tau=0.25)
QR_75 <- rq(IntStr ~ TL, tau=0.75)
anova(QR_25, QR_75)  # test difference btw quantiles 25 & 75

# Scatter plot with quantile regression lines

q <- c(0.25, 0.75)

qr_IS_TL <- ggplot(spp_attr_all, aes(x = TLu, y = AllStrength_mean)) +
  # geom_errorbar(aes(ymin = AllStrength_mean-AllStrength_Q1, ymax = AllStrength_mean+AllStrength_Q3), width = .1, position = pd) +
  geom_point() + 
  geom_quantile(quantiles = q, size = 2, alpha = 0.5, aes(colour = as.factor(..quantile..))) +
  scale_y_log10() +
  labs(x = "Trophic level", y = "mean Interaction Strength", color = "Quantiles") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
qr_IS_TL


## Degree ----




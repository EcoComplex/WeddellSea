#
## Fit distribution of interaction strength using Maximum Likelihood
## 'univariateML' package
#


## Load packages

packages <- c("univariateML", "fitdistrplus", "gamlss", "ggplot2", "Rmisc", "dplyr")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


## Load data

load("Results/network_&_spp_attr.rda")


## Distribution of interaction strength (qRC)

# Plot
ggplot(wedd_int, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  scale_y_log10()


# Distribution fit

x <- wedd_int$qRC
aic.result <- c(AIC(mlunif(x), mlexp(x), mlpower(x), mllnorm(x), mlnorm(x), mlgamma(x)))
IS_fit <- bind_cols(aic.result) %>% 
  mutate(Model = c("Uniform", "Exponential", "Power-law", "log-Normal", "Normal", "Gamma"),
         deltaAIC = AIC - min(AIC)) %>% 
  arrange(deltaAIC) %>% 
  dplyr::select(Model, df, AIC, deltaAIC)

IS_fit


#
# Power law with exponential cutoff 
#
source("R/evaluate_distr.r")

pdist <- evaluate_distr(wedd_int$qRC,est_xmin=FALSE,returnOBJ=TRUE)
pdist1 <- evaluate_distr(wedd_int$qRC,est_xmin=FALSE,returnOBJ=FALSE)

dist_manejo %>% filter(min(AICc)==AICc) %>% mutate(moda=exp(expo))


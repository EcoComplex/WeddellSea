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


# Function to fit distribution & test it
aic.bic <- function(x){
  
  # Akaike's Information Criterion
  aic.result <- c(AIC(mlunif(x), mlexp(x), mlpower(x), mllnorm(x), mlnorm(x), mlgamma(x)))
  
  # Schwarz's Bayesian criterion
  bic.result <- c(BIC(mlunif(x), mlexp(x), mlpower(x), mllnorm(x), mlnorm(x), mlgamma(x)))
  
  aic.result <- aic.result$AIC
  bic.result <- bic.result$BIC
  names(aic.result) <- c("uniform", "exp", "power", "lnorm", "norm", "gamma")
  names(bic.result) <- c("uniform", "exp", "power", "lnorm", "norm", "gamma")
  
  result <- bind_rows(aic.result, bic.result) %>% 
    mutate(Statistic = c("AIC", "BIC"), 
           BestFit = c(names(aic.result)[which.min(aic.result)], names(bic.result)[which.min(bic.result)])) %>% 
    dplyr::select(Statistic, everything())
  
  return(result)
}

IS_fit <- aic.bic(wedd_int$qRC)
IS_fit

descdist(wedd_int$qRC, discrete = TRUE)

#
# Power law with exponential cutoff 
#
source("R/evaluate_distr.r")

pdist <- evaluate_distr(wedd_int$qRC,est_xmin=FALSE,returnOBJ=TRUE)
pdist1 <- evaluate_distr(wedd_int$qRC,est_xmin=FALSE,returnOBJ=FALSE)

dist_manejo %>% filter(min(AICc)==AICc) %>% mutate(moda=exp(expo))


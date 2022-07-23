#
## Fit distribution of interaction strength using Maximum Likelihood
## 'univariateML' package
#


## Load packages

packages <- c("univariateML", "fitdistrplus", "gamlss", "ggplot2", "Rmisc")
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
  
  aic.result <- c(AIC(mlunif(x), mlexp(x), mlpower(x), mllnorm(x), mlnorm(x)))
  
  bic.result <- c(BIC(mlunif(x), mlexp(x), mlpower(x), mllnorm(x), mlnorm(x)))
  
  aic.result <- aic.result$AIC
  bic.result <- bic.result$BIC
  names(aic.result) <- c("uniform", "exp", "power", "lnorm", "norm")
  names(bic.result) <- c("uniform", "exp", "power", "lnorm", "norm")
  
  result <- c(names(aic.result)[which.min(aic.result)], names(bic.result)[which.min(bic.result)])
  
  return(result)
}

aic.bic(wedd_int$qRC)
descdist(wedd_int$qRC, discrete = TRUE)


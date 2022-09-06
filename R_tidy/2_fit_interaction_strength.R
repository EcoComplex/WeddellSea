#
## Fit distribution of interaction strength estimation
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first run 'R/1_calc_interaction_strength.R'


## Load packages ----

packages <- c("univariateML", "fitdistrplus", "gamlss", "ggplot2", "Rmisc", "dplyr")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


## Load interaction strength estimation ----

load("Results/interaction_estimation.rda")


## Distribution of interaction strength ----

# Plot
ggplot(wedd_int_pd, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  scale_y_log10() +
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))


# Fit distribution using Maximum Likelihood
x <- wedd_int_pd$qRC
aic.result <- c(AIC(mlunif(x), mlexp(x), mlpower(x), mllnorm(x), mlnorm(x), mlgamma(x)))
IS_fit <- bind_cols(aic.result) %>% 
  mutate(Model = c("Uniform", "Exponential", "Power-law", "log-Normal", "Normal", "Gamma"),
         deltaAIC = AIC - min(AIC)) %>% 
  arrange(deltaAIC) %>% 
  dplyr::select(Model, df, AIC, deltaAIC)
IS_fit

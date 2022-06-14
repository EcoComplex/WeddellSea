# 
# Build a GAM model to predict interactionStrength from topological measures
#

require(tidyverse)
spp_attr <- read_csv("Data/spp_all_attr.csv")



## GAM prediction of interaction strength with TL and degree ----
#
require(mgcv)
require(gratia)

# Full model 

fit_full <- gam(log(TotalStrength) ~ te(TLu, Degree, meanTrophicSimil, k=c(6,6,6)), data = spp_attr)
summary(fit_full)
gam.check(fit_full)
draw(fit_full)
appraise(fit_full) 


# Use weighted TL
fitw <- gam(log(TotalStrength) ~ te(TLw, Degree, k=c(8,8)), data = spp_attr)
plot(fitw, rug=F, pers=T, theta=45, main="Strength")
draw(fitw, residuals=T) 
appraise(fitw) 
gam.check(fitw)
summary(fitw)  # percentage of deviance explained


fit_dt <- gam(log(TotalStrength) ~ te( Degree, meanTrophicSimil, k=c(6,6)), data = spp_attr)
summary(fit_dt)  # percentage of deviance explained

fit_tt <- gam(log(TotalStrength) ~ te( TLu, meanTrophicSimil, k=c(6,6)), data = spp_attr)
summary(fit_tt)  # percentage of deviance explained

# Use unweighted (topological) TL
fitu <- gam(log(TotalStrength) ~ te(TLu, Degree, k=c(8,8)), data = spp_attr)
plot(fitu, rug=F, scheme=T, theta=45, main="Strength")
draw(fitu, residuals=T) 
appraise(fitu) 
gam.check(fitu)
summary(fitu)

# Compare GAM models (using weighted & unweighted TL)

delta_aic <- tibble(AIC(fit_full,fit_dt,fit_tt,fitu,fitw))
delta_aic$model <- rownames(AIC(fit_full,fit_dt,fit_tt,fitu,fitw))

delta_aic %>% mutate(delta_AIC = AIC - min(AIC) ) %>% arrange(delta_AIC)





# Plot GAM model and observed data
# Taken from https://stackoverflow.com/questions/55047365/r-plot-gam-3d-surface-to-show-also-actual-response-values
# 2D plot

fitu_plot <- gam(TotalStrength ~ te(TLu, Degree, k=c(6,6)), data = spp_attr)  # not work with 'family=tw'
summary(fitu_plot)  # deviance % explained by the model
appraise(fitu_plot) 

# Now expand it to a grid so that persp will work
steps <- 30
Degree <- with(spp_attr, seq(min(Degree), max(Degree), length = steps) )
TLu <-  with(spp_attr, seq(min(TLu), max(TLu), length = steps) )
newdat <- expand.grid(Degree = Degree, TLu = TLu)
ilink <- family(fitu_plot)$linkinv 

TotalStrength <- matrix(predict(fitu_plot, newdat), steps, steps)
TotalStrength <- ilink(TotalStrength)

# Now plot it
p <- persp(Degree, TLu, TotalStrength, theta = 45, col = alpha("yellow", 0.2), 
           xlim = range(Degree), ylim = range(TLu), zlim = range(TotalStrength),
           xlab = "Degree", ylab = "Trophic level", zlab = "Interaction strength")

# To add the points, you need the same 3d transformation
obs <- with(spp_attr, trans3d(Degree, TLu, TotalStrength, p))
pred <- with(spp_attr, trans3d(Degree, TLu, fitted(fitu_plot), p))
dif <- obs[["y"]] - pred[["y"]]  # observed - predicted strength
points(obs, col = ifelse(dif < 0, alpha("blue",0.4), alpha("red",1)), pch = 16)
# Add segments to show where the points are in 3d
segments(obs$x, obs$y, pred$x, pred$y)

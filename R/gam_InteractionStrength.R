# 
## Build a GAM model to describe the dependence of interaction strength
## using trophic level, degree and trophic similarity
#


# Load packages ----

require(tidyverse)
require(mgcv)
require(gratia)


# Load data ----

load("Results/network_&_spp_attr.rda")


# General Additive Model for mean interaction strength ----

## Full model ----

# Trophic level, degree & mean trophic similarity
fit_full <- gam(AllStrength_mean ~ te(TLu, Degree, meanTrophicSimil, k=c(6,6,6)), data = spp_attr_all)
summary(fit_full)
gam.check(fit_full)
draw(fit_full)
appraise(fit_full) 

# Full model no interaction
fit_fnint <- gam(AllStrength_mean ~ s(TLu) + s(Degree) + s(meanTrophicSimil), data = spp_attr_all)
summary(fit_fnint)
gam.check(fit_fnint)
draw(fit_fnint)
appraise(fit_fnint) 

## TL & Degree ----

fit_tldeg <- gam(AllStrength_mean ~ te(TLu, Degree, k=c(6,6)), data = spp_attr_all)
summary(fit_tldeg)
gam.check(fit_tldeg)
draw(fit_tldeg)
appraise(fit_tldeg) 

# TL & Degree no interaction
fit_tldeg_nint <- gam(AllStrength_mean ~ s(TLu) + s(Degree), data = spp_attr_all)
summary(fit_tldeg_nint)
gam.check(fit_tldeg_nint)
draw(fit_tldeg_nint)
appraise(fit_tldeg_nint) 

## TL & TS ----

fit_tlts <- gam(AllStrength_mean ~ te(TLu, meanTrophicSimil, k=c(6,6)), data = spp_attr_all)
summary(fit_tlts)
gam.check(fit_tlts)
draw(fit_tlts)
appraise(fit_tlts) 

# TL & TS no interaction
fit_tlts_nint <- gam(AllStrength_mean ~ s(TLu) + s(meanTrophicSimil), data = spp_attr_all)
summary(fit_tlts_nint)
gam.check(fit_tlts_nint)
draw(fit_tlts_nint)
appraise(fit_tlts_nint)

## Degree & TS ----

fit_degts <- gam(AllStrength_mean ~ te(Degree, meanTrophicSimil, k=c(6,6)), data = spp_attr_all)
summary(fit_degts)
gam.check(fit_degts)
draw(fit_degts)
appraise(fit_degts) 

# TL & TS no interaction
fit_degts_nint <- gam(AllStrength_mean ~ s(Degree) + s(meanTrophicSimil), data = spp_attr_all)
summary(fit_degts_nint)
gam.check(fit_degts_nint)
draw(fit_degts_nint)
appraise(fit_degts_nint)


# Compare GAMs ----

delta_aic <- tibble(AIC(fit_full, fit_fnint, fit_tldeg, fit_tldeg_nint,
                        fit_tlts, fit_tlts_nint, fit_degts, fit_degts_nint))
delta_aic$model <- rownames(AIC(fit_full, fit_fnint, fit_tldeg, fit_tldeg_nint,
                                fit_tlts, fit_tlts_nint, fit_degts, fit_degts_nint))

delta_aic %>% mutate(delta_AIC = AIC - min(AIC) ) %>% arrange(delta_AIC)


# Use weighted TL
fitw <- gam(TotalStrength ~ te(TLw, Degree, k=c(8,8)), data = spp_attr)
plot(fitw, rug=F, pers=T, theta=45, main="Strength")
draw(fitw, residuals=T) 
appraise(fitw) 
gam.check(fitw)
summary(fitw)  # percentage of deviance explained


fit_dt <- gam(TotalStrength ~ te( Degree, meanTrophicSimil, k=c(6,6)), data = spp_attr)
summary(fit_dt)  # percentage of deviance explained
appraise(fit_dt) 

fit_tt <- gam(TotalStrength ~ te( TLu, meanTrophicSimil, k=c(6,6)), data = spp_attr)
summary(fit_tt)  # percentage of deviance explained
appraise(fit_dt) 

# Use unweighted (topological) TL
fitu <- gam(TotalStrength ~ te(TLu, Degree, k=c(8,8)), data = spp_attr)
plot(fitu, rug=F, scheme=T, theta=45, main="Strength")
draw(fitu, residuals=T) 
appraise(fitu) 
gam.check(fitu)
summary(fitu)







# Plot GAM model and observed data
# Taken from https://stackoverflow.com/questions/55047365/r-plot-gam-3d-surface-to-show-also-actual-response-values
# 2D plot
fitu_plot <- gam((TotalStrength) ~ te(TLu, Degree, k=c(6,6)), data = spp_attr,family=tw)

# fitu_plot <- gam(TotalStrength ~ te(TLu, Degree,), data = spp_attr)  # not work with 'family=tw'
summary(fitu_plot)  # deviance % explained by the model
appraise(fitu_plot) 

# Now expand it to a grid so that persp will work
steps <- 30
Degree <- with(spp_attr, seq(min(Degree), max(Degree), length = steps) )
TLu <-  with(spp_attr, seq(min(TLu), max(TLu), length = steps) )
#meanTrophicSimil <- with(spp_attr, seq(min(meanTrophicSimil), max(meanTrophicSimil), length = 6) )
newdat <- expand.grid(Degree = Degree, TLu = TLu)
ilink <- family(fitu_plot)$linkinv 

TotalStrength <- matrix(predict(fitu_plot, newdat), steps, steps)
TotalStrength <- ilink(TotalStrength)

# Now plot it
p <- persp(Degree, TLu, TotalStrength, theta = 10, col = alpha("yellow", 0.2), 
           xlim = range(Degree), ylim = range(TLu), zlim = range(TotalStrength),
           xlab = "Degree", ylab = "Trophic level", zlab = "Interaction strength")

# To add the points, you need the same 3d transformation
obs <- with(spp_attr , trans3d(Degree, TLu, (TotalStrength), p))
pred <- with(spp_attr, trans3d(Degree, TLu, fitted(fitu_plot), p))
dif <- obs[["y"]] - pred[["y"]]  # observed - predicted strength
points(obs, col = ifelse(dif < 0, alpha("blue",0.4), alpha("red",1)), pch = 16)
# Add segments to show where the points are in 3d
segments(obs$x, obs$y, pred$x, pred$y)

require(plotly)
plot_ly() %>% 
  add_surface(x = ~Degree, y = ~TLu, z = ~TotalStrength) %>% 
  add_markers(x = ~Degree, y=~TLu, z=~(TotalStrength),data=spp_attr,size=1,hoverinfo = "text",
              text = ~paste(TrophicSpecies, '<br>Degree:', Degree, '<br>Trophic level:', TLu, 
                            '<br>Interaction strength:', round(TotalStrength,4))) %>% 
  layout(scene = list(
          yaxis = list(autorange="reversed"),
          xaxis = list(autorange="reversed"),
          zaxis = list(range = list(0, 0.1))
          )
  )
  

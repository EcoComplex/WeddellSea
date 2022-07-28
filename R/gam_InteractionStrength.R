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
de_full <- summary(fit_full)$dev.expl
gam.check(fit_full)
draw(fit_full)
appraise(fit_full) 

# Full model no interaction
fit_fnint <- gam(AllStrength_mean ~ s(TLu) + s(Degree) + s(meanTrophicSimil), data = spp_attr_all)
summary(fit_fnint)
de_fnint <- summary(fit_fnint)$dev.expl
gam.check(fit_fnint)
draw(fit_fnint)
appraise(fit_fnint) 

## TL & Degree ----

fit_tldeg <- gam(AllStrength_mean ~ te(TLu, Degree, k=c(6,6)), data = spp_attr_all)
summary(fit_tldeg)
de_tldeg <- summary(fit_tldeg)$dev.expl
gam.check(fit_tldeg)
draw(fit_tldeg)
appraise(fit_tldeg)
plot(fit_tldeg, rug=F, scheme=T, theta=45, 
     main="Interaction Strength", xlab="Trophic level", ylab="Degree")

# TL & Degree no interaction
fit_tldeg_nint <- gam(AllStrength_mean ~ s(TLu) + s(Degree), data = spp_attr_all)
summary(fit_tldeg_nint)
de_tldeg_nint <- summary(fit_tldeg_nint)$dev.expl
gam.check(fit_tldeg_nint)
draw(fit_tldeg_nint)
appraise(fit_tldeg_nint) 

## TL & TS ----

fit_tlts <- gam(AllStrength_mean ~ te(TLu, meanTrophicSimil, k=c(6,6)), data = spp_attr_all)
summary(fit_tlts)
de_tlts <- summary(fit_tlts)$dev.expl
gam.check(fit_tlts)
draw(fit_tlts)
appraise(fit_tlts) 

# TL & TS no interaction
fit_tlts_nint <- gam(AllStrength_mean ~ s(TLu) + s(meanTrophicSimil), data = spp_attr_all)
summary(fit_tlts_nint)
de_tlts_nint <- summary(fit_tlts_nint)$dev.expl
gam.check(fit_tlts_nint)
draw(fit_tlts_nint)
appraise(fit_tlts_nint)

## Degree & TS ----

fit_degts <- gam(AllStrength_mean ~ te(Degree, meanTrophicSimil, k=c(6,6)), data = spp_attr_all)
summary(fit_degts)
de_degts <- summary(fit_degts)$dev.expl
gam.check(fit_degts)
draw(fit_degts)
appraise(fit_degts) 

# TL & TS no interaction
fit_degts_nint <- gam(AllStrength_mean ~ s(Degree) + s(meanTrophicSimil), data = spp_attr_all)
summary(fit_degts_nint)
de_degts_nint <- summary(fit_degts_nint)$dev.expl
gam.check(fit_degts_nint)
draw(fit_degts_nint)
appraise(fit_degts_nint)


# Compare GAMs ----

delta_aic <- tibble(AIC(fit_full, fit_fnint, fit_tldeg, fit_tldeg_nint,
                        fit_tlts, fit_tlts_nint, fit_degts, fit_degts_nint))
delta_aic$model <- rownames(AIC(fit_full, fit_fnint, fit_tldeg, fit_tldeg_nint,
                                fit_tlts, fit_tlts_nint, fit_degts, fit_degts_nint))
delta_aic$dev.expl <- c(de_full, de_fnint, de_tldeg, de_tldeg_nint, de_tlts,
                        de_tlts_nint, de_degts, de_degts_nint)
delta_aic %>% mutate(delta_AIC = AIC - min(AIC) ) %>% 
  arrange(delta_AIC) %>% dplyr::select(model, dev.expl, AIC, df, delta_AIC)


# Plot GAM ----

str(spp_attr_all)
v <- ggplot(spp_attr_all, aes(x = Degree, y = TLu, z = AllStrength_mean))
v + geom_contour_filled()
v + geom_raster(aes(fill = AllStrength_mean)) +
  geom_contour(colour = "white")

plot.gam(fit_full, pages=1, residuals=TRUE, all.terms=TRUE, shade=TRUE, shade.col=2)
plot.gam(fit_full, pages=1, seWithMean=TRUE)

library(mgcViz)
b <- getViz(fit_full)
plotRGL(sm(b, 1), fix = c("z" = 0), residuals = TRUE)

str(faithfuld)
g <- ggplot(faithfuld, aes(waiting, eruptions, z = density))
g + geom_contour_filled()





# Plot GAM model and observed data
# Taken from https://stackoverflow.com/questions/55047365/r-plot-gam-3d-surface-to-show-also-actual-response-values
# 2D plot
fitu_plot <- gam((AllStrength_mean) ~ te(TLu, Degree, k=c(6,6)), data = spp_attr_all)  # not work with 'family=tw'
summary(fitu_plot)  # deviance % explained by the model
appraise(fitu_plot) 

# Now expand it to a grid so that persp will work
steps <- 30
Degree <- with(spp_attr_all, seq(min(Degree), max(Degree), length = steps) )
TLu <-  with(spp_attr_all, seq(min(TLu), max(TLu), length = steps) )
#meanTrophicSimil <- with(spp_attr_all, seq(min(meanTrophicSimil), max(meanTrophicSimil), length = 6))
newdat <- expand.grid(Degree = Degree, TLu = TLu)
ilink <- family(fitu_plot)$linkinv 

AllStrength_mean <- matrix(predict(fitu_plot, newdat), steps, steps)
AllStrength_mean <- ilink(AllStrength_mean)

# Now plot it
p <- persp(Degree, TLu, AllStrength_mean, theta = 10, col = alpha("yellow", 0.2), 
           xlim = range(Degree), ylim = range(TLu), zlim = range(AllStrength_mean),
           xlab = "Degree", ylab = "Trophic level", zlab = "mean Interaction strength")

# To add the points, you need the same 3d transformation
obs <- with(spp_attr_all, trans3d(Degree, TLu, (AllStrength_mean), p))
pred <- with(spp_attr_all, trans3d(Degree, TLu, fitted(fitu_plot), p))
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
  

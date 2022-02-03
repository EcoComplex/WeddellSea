## Estimate interaction strength
# Author: TIM
# Date: January 2022


## Load packages ----

pkg_list <- c("dplyr", "readr", "ggplot2", "ggpmisc", "IDPmisc")

lapply(pkg_list, FUN = function(pkg_list) {
  do.call("require", list(pkg_list)) 
})


## Load data ----

# Load total (data_total) & curated data (weddell_df) from 1_datawrangling.R
load("../Data/Weddell_database.rda")
load("../Data/g_unweighted_att.rda")


## Completing missing interaction dimensionality ----

# Following Pawar et al. 2012 criteria:
# benthic resource (R) and benthic consumer (C) = 2D
# benthic resource and pelagic consumer = 2D
# pelagic resource and pelagic consumer = 3D
# pelagic resource and benthic consumer = 3D

# Then, looking at res.movement.type and con.movement.type we can get dimensionality:
# 'sessile' R and 'sessile' C = 2D
# 'sessile' R and 'walking' C = 2D
# 'walking' R and 'walking' C = 2D
# 'walking' R and 'sessile' C = 2D
# 'swimming' R and any con.movement.type C = 3D

weddell_dim_fill <- weddell_df %>% 
  select(res.taxonomy, res.movement.type, res.mass.mean.g.,
         con.taxonomy, con.movement.type, con.mass.mean.g., interaction.dimensionality) %>% 
  mutate (complete.dim = case_when(res.movement.type == "sessile" ~ "2D",
                               res.movement.type == "walking" ~ "2D",
                               res.taxonomy == "Sediment" ~ "2D",
                               res.movement.type == "swimming" ~ "3D"))

weddell_dim_fill <- weddell_dim_fill %>%
  select(con.taxonomy, res.taxonomy, con.mass.mean.g., res.mass.mean.g., complete.dim) %>% 
  rename(interaction.dim = complete.dim)
write.table(weddell_dim_fill, file = "../Data/Wedd_mass_complete.dat")


## Estimation of search rate "alfa" (following Pawar et al. 2012) ----

# alfa_2D = alfa2D * mC^0.68, where alfa2D = -3.08
# alfa_3D = alfa3D * mC^1.05, where alfa3D = 1.77
# derivative alfa_2D = alpha2D * ln(mC)*mC^0.68
# derivative alfa_3D = alpha3D * ln(mC)*mC^1.05

# Estimate derivative alfa_2D & alfa_3D
weddell_sr <- weddell_df %>% 
  mutate(m_R.kg = res.mass.mean.g.*1e-3, m_C.kg = con.mass.mean.g.*1e-3) %>%  # g to kg
  mutate(mass.ratio = m_R.kg/m_C.kg) %>% 
  mutate(alfa_2D = -3.08*m_C.kg^0.68, alfa_3D = 1.77*m_C.kg^1.05) %>%  # search rate
  mutate(der_alfa_2D = -3.08*(log(m_C.kg)*m_C.kg^0.68), der_alfa_3D = -1.77*(log(m_C.kg)*m_C.kg^1.05))  # derivative search rate

ggplot(weddell_sr, aes(x = log10(m_C.kg), y = log10(m_R.kg))) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA)

# Fit to ordinary least squares regression (OLS) & plot m_C vs der_alfa
weddell_OLS <- weddell_sr %>% 
  mutate(log_der_alfa_2D = log10(der_alfa_2D), log_der_alfa_3D = log10(der_alfa_3D), 
         log_m_C.kg = log10(m_C.kg))
weddell_OLS <- NaRV.omit(weddell_OLS)  # omit NA, NaN, Inf

# Interaction dimensionality 2D
OLS_model_2D <- lm(log_der_alfa_2D ~ log_m_C.kg, data = weddell_OLS)
summary(OLS_model_2D)

plot_OLS_2D <- ggplot(data = weddell_OLS, aes(x = log_m_C.kg, y = log_der_alfa_2D)) +
  geom_point(colour = "red") +
  geom_smooth(method = "lm", se=FALSE, color="black", linetype = "dashed", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  xlab("log(m_C.kg)") + ylab("log(der_alfa_2D)")
plot_OLS_2D

# Interaction dimensionality 3D
OLS_model_3D <- lm(log_der_alfa_3D ~ log_m_C.kg, data = weddell_OLS)
summary(OLS_model_3D)

plot_OLS_3D <- ggplot(data = weddell_OLS, aes(x = log_m_C.kg, y = log_der_alfa_3D)) +
  geom_point(colour = "blue") +
  geom_smooth(method = "lm", se=FALSE, color="black", linetype = "dashed", formula = y ~ x) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  xlab("log(m_C.kg)") + ylab("log(der_alfa_3D)")
plot_OLS_3D


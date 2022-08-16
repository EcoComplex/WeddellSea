#
## Estimate interaction strength following Pawar et al. (2012)
## 'multiweb' package
#


# Load packages ----

packages <- c("dplyr", "readr", "multiweb", "ggplot2")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


# Load food web database ----

# Read and check data available on GATEWAy about Weddell Sea
# Read the original database
ga <- readr::read_csv("Data/283_2_FoodWebDataBase_2018_12_10.csv", col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
names(ga)
wedd_df <-  ga %>%
  filter(grepl('Weddell', foodweb.name))

# Rename and select fields
wedd_df <- wedd_df %>% 
  rename(con_mass_mean = con.mass.mean.g., res_mass_mean = res.mass.mean.g.,
         con_taxonomy = con.taxonomy, res_taxonomy = res.taxonomy, interaction_dim = interaction.dimensionality) %>% 
  dplyr::select(con_taxonomy, res_taxonomy, con_mass_mean,res_mass_mean, interaction_dim)


# Completing missing interaction dimensionality ----

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

weddell_dim_fill <- wedd_df %>% 
  select(res.taxonomy, res.movement.type, res.mass.mean.g.,
         con.taxonomy, con.movement.type, con.mass.mean.g., interaction.dimensionality) %>% 
  mutate (complete.dim = case_when(res.movement.type == "sessile" ~ "2D",
                               res.movement.type == "walking" ~ "2D",
                               res.taxonomy == "Sediment" ~ "2D",
                               res.movement.type == "swimming" ~ "3D"))
weddell_dim_fill <- weddell_dim_fill %>%
  select(con.taxonomy, res.taxonomy, con.mass.mean.g., res.mass.mean.g., complete.dim) %>% 
  rename(interaction.dim = complete.dim)

write.table(weddell_dim_fill, file = "Data/Wedd_mass_complete.dat")


# Estimation of interaction strength ----

# Read the updated database
# We completed missing interaction dimensionality following Pawar et al. (2012)

wedd_df_comp <- read_delim("Data/Wedd_mass_complete.dat", delim = " ", col_types="dccddc") %>% 
  mutate(res_mass_mean = res_mass_mean*1000, con_mass_mean = con_mass_mean*1000)

# Replace 'phytodetritus' and 'sediment' body masses (-999) with smallest phytoplankton
# Smallest phyto body mass = 1.53e-11
wedd_df_pd <- wedd_df_comp %>% 
  mutate(res_mass_mean = replace(res_mass_mean, res_mass_mean == -999000, 1.53e-11))

# Estimate interaction strength
wedd_int_pd <- multiweb::calc_interaction_intensity(wedd_df_pd, res_mass_mean, con_mass_mean, interaction_dim)

# Explore distribution of interaction intensity (qRC)
ggplot(wedd_int_pd, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  scale_y_log10() +
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))


# Dependence of IS estimation ----

# On mR/mC
plot <- wedd_int_pd %>% 
  mutate(mR_mC = res_mass_mean/con_mass_mean) %>% 
  ggplot(., aes(x = xR, y = qRC)) +
  geom_point()
plot

# On mR & mC 3D
require(plotly)
plot_3D <- wedd_int_pd %>% 
  mutate(mR_mC = res_mass_mean/con_mass_mean) %>%
  filter(interaction_dim == "3D") %>% 
  plot_ly(data=., x=~res_mass_mean, y=~con_mass_mean, 
                   z=~log(qRC), type="scatter3d", mode="markers", 
                   color="black", marker = list(size = 5), hoverinfo = "text",
                   text = ~paste('<br>Res_mass:', res_mass_mean, '<br>Con_mass:', con_mass_mean, 
                                 '<br>Interaction strength:', round(qRC,5))) %>%  
  layout(scene = list(xaxis=list(title = "Res_mass"),
                      yaxis=list(title = "Con_mass"),
                      zaxis=list(title = "log Interaction strength", color = "red")))
plot_3D

# On mR & mC 2D
plot_2D <- wedd_int_pd %>% 
  mutate(mR_mC = res_mass_mean/con_mass_mean) %>%
  filter(interaction_dim == "2D") %>% 
  plot_ly(data=., x=~res_mass_mean, y=~con_mass_mean, 
          z=~log(qRC), type="scatter3d", mode="markers", 
          color="black", marker = list(size = 5), hoverinfo = "text",
          text = ~paste('<br>Res_mass:', res_mass_mean, '<br>Con_mass:', con_mass_mean, 
                        '<br>Interaction strength:', round(qRC,5))) %>%  
  layout(scene = list(xaxis=list(title = "Res_mass"),
                      yaxis=list(title = "Con_mass"),
                      zaxis=list(title = "log Interaction strength", color = "red")))
plot_2D


# Save results ----

save(wedd_df_pd, wedd_int_pd,
     file = "Results/interaction_estimation.rda")

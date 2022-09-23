#
## Estimate interaction strength using 'multiweb' package 0.5.20
## Authors: Leonardo Saravia & Tomás Ignacio Marina
## September 2022
#

## To run the following code first clone the GitHub repo https://github.com/EcoComplex/WeddellSea
## git clone https://github.com/EcoComplex/WeddellSea.git


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

# Read and check data available on GATEWAy version 1.0
# https://doi.org/10.25829/IDIV.283-3-756
# Filter Weddell Sea food web database
ga <- readr::read_csv("Data/283_2_FoodWebDataBase_2018_12_10.csv", col_types = "icccccccccccccddddddccccccccddddddcddccddciicc")
names(ga)
wedd_df <-  ga %>%
  filter(grepl('Weddell', foodweb.name))


# Complete missing interaction dimensionality ----

# Following Pawar et al. 2012 (https://doi.org/10.1038/nature11131) criteria:
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
  rename(interaction.dimensionality = complete.dim)

#write_csv(weddell_dim_fill, file = "Data/Wedd_int_complete.csv")


# Estimation of interaction strength ----

# Read the updated database
# Convert g to kg (to follow Pawar et al. 2012 relationships)
wedd_df_comp <- read_csv("Data/Wedd_int_complete.csv") %>% 
  mutate(res.mass.mean.kg. = res.mass.mean.g.*10e-3, con.mass.mean.kg. = con.mass.mean.g.*10e-3)

# Replace 'phytodetritus' and 'sediment' body masses (-999) with that of smallest phytoplankton
# Smallest phyto 'Fragilariopsis cylindrus' mean body mass (1.53e-14 g) * 10e-3 = 1.5e-16 kg
wedd_df_pd <- wedd_df_comp %>% 
  mutate(res.mass.mean.kg. = replace(res.mass.mean.kg., res.mass.mean.kg. == -999000, 1.53e-16))

# Estimate interaction strength
wedd_int_pd <- multiweb::calc_interaction_intensity(wedd_df_pd, res.mass.mean.kg., con.mass.mean.kg., interaction.dimensionality)

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


# Save results ----

save(wedd_df_pd, wedd_int_pd,
     file = "Results/interaction_estimation.rda")

#
## Figures showing relationship btw Interaction Strength 'IS' & spp attributes
#


# Load packages ----

library(ggplot2)
library(tidyverse)
library(Rmisc)
library(ggtext)
library(plotly)
source("R/network_fun.r")

# Load spp data ----

load("Results/network_&_spp_attr.rda")


# Explore distribution of interaction intensity (qRC) ----

ggplot(wedd_int, aes(log(qRC))) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18))


# Interaction Strength & spp attr. ---- 

# Cumulative strength & proportion of spp 

spp_attr_all <- spp_attr_all %>% 
  mutate(cumsum_str = order_by(-AllStrength_sum, cumsum(AllStrength_sum)),
         prop_str = cumsum_str/(sum(AllStrength_sum))) %>% 
  mutate(rank_spp = dense_rank(desc(AllStrength_sum)),
         prop_spp = rank_spp/nrow(spp_attr_all))

## By Trophic Level ----

summary(spp_attr_all$TLu)
spp_nt_12 <- length(spp_attr_all$TLu[spp_attr_all$TLu < 2])
spp_nt_23 <- length(spp_attr_all$TLu[spp_attr_all$TLu >= 2 & spp_attr_all$TLu < 3])
spp_nt_34 <- length(spp_attr_all$TLu[spp_attr_all$TLu >= 3 & spp_attr_all$TLu < 4])
spp_nt_45 <- length(spp_attr_all$TLu[spp_attr_all$TLu >= 4 & spp_attr_all$TLu < 5])
spp_nt_56 <- length(spp_attr_all$TLu[spp_attr_all$TLu < 6])

pd <- position_dodge(0.5) # move them .05 to the left and right

IS_TL <- spp_attr_all %>% 
  # mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>% 
  ggplot(aes(x = reorder(TrophicSpecies, TLu), y = AllStrength_mean)) + 
  geom_errorbar(aes(ymin = AllStrength_mean-AllStrength_Q1, ymax = AllStrength_mean+AllStrength_Q3), width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +  # aes(color = Color)
  scale_y_log10() +
  # scale_color_identity() +
  geom_vline(xintercept = c(spp_nt_12+1, spp_nt_12+spp_nt_23+1, spp_nt_12+spp_nt_23+spp_nt_34+1, 
                            spp_nt_12+spp_nt_23+spp_nt_34+spp_nt_45+1, spp_nt_56+1), 
             linetype = "longdash", colour = "blue") +
  labs(x = "Trophic species (ascending TL)", y = "Interaction strength (log scale)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15))

IS_TL + 
  annotate(x = spp_nt_12+1, y = +Inf, label = "Trophic Level 2", geom = "label", vjust = 2) +
  annotate(x = spp_nt_12+spp_nt_23+1, y = +Inf, label = "Trophic Level 3", geom = "label", vjust = 2) +
  annotate(x = spp_nt_12+spp_nt_23+spp_nt_34+1, y = +Inf, label = "Trophic Level 4", geom = "label", vjust = 2) # +
  # annotate(x = spp_nt_12+spp_nt_23+spp_nt_34+20, y = +Inf, label = "Trophic Level 5", geom = "label", vjust = 8)


## By Degree ----

IS_Deg <- spp_attr_all %>% 
  mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>% 
  ggplot(aes(x = reorder(TrophicSpecies, Degree), y = AllStrength_mean)) + 
  geom_errorbar(aes(ymin = AllStrength_mean-AllStrength_Q1, ymax = AllStrength_mean+AllStrength_Q3), width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, aes(color = Color)) +
  scale_y_log10() +
  scale_color_identity() +
  labs(x = "Trophic species (ascending degree)", y = "Interaction strength (log scale)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15))
IS_Deg

ggplot(spp_attr_all, aes(x = Degree, y = AllStrength_mean)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Degree (log scale)", y = "Interaction strength (log scale)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

## By Trophic Similarity ----

IS_TS <- spp_attr_all %>% 
  mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>% 
  ggplot(aes(x = reorder(TrophicSpecies, -meanTrophicSimil), y = AllStrength_mean)) +
  geom_errorbar(aes(ymin = AllStrength_mean-AllStrength_Q1, ymax = AllStrength_mean+AllStrength_Q3), width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, aes(color = Color)) +
  scale_y_log10() +
  scale_color_identity() +
  labs(x = "Trophic species (descending Trophic simil.)", y = "Interaction strength (log scale)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15))
IS_TS


# Interactive plot (IS, Degree, TL) ----

plot_3d <- plot_ly(data=spp_attr_all, x=~Degree, y=~TLu, 
                   z=~AllStrength_mean, type="scatter3d", mode="markers", 
                   color="black", marker = list(size = 5), hoverinfo = "text",
                   text = ~paste(TrophicSpecies, '<br>Degree:', Degree, '<br>Trophic level:', TLu, 
                                 '<br>Mean IS:', round(AllStrength_mean,4))) %>%  
  layout(scene = list(xaxis=list(title = "Degree"),
                      yaxis=list(title = "Trophic level"),
                      zaxis=list(title = "Interaction strength", color = "red")))
plot_3d


# Trophic Level by T. Similarity ----

TL_TS <- spp_attr_all %>% 
  mutate(Color = ifelse(prop_str < 0.8, "red", "black")) %>% 
  ggplot(aes(x = meanTrophicSimil, y = TLu)) +
  geom_point(aes(color = Color)) +
  scale_color_identity() +
  labs(x = "Trophic similarity", y = "Trophic level") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
TL_TS

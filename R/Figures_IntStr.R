#
## Figures showing relationship btw Interaction Strength 'IS' & sp properties
#

# Load packages
library(ggplot2)
library(tidyverse)
library(Rmisc)
library(ggtext)
library(plotly)
source("R/network_fun.r")

# Load spp data
load("Data/network_&_spp_attr.rda")


# Explore distribution of interaction intensity (qRC)
ggplot(wedd_int, aes(qRC)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "white") + 
  labs(x = "Interaction strength", y = "Frequency (log scale)") +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) +
  scale_y_log10()


# Summary of IS
sum_IS <- summarySE(all_int, measurevar="IS", groupvars="TrophicSpecies")
sum_IS

# Join spp attributes to 'sum_IS'
all_data <- spp_attr_all %>% 
  left_join(sum_IS)


## Plot IS by Trophic Level

summary(all_data$TLu)
spp_nt_12 <- length(all_data$TLu[all_data$TLu < 2])
spp_nt_23 <- length(all_data$TLu[all_data$TLu >= 2 & all_data$TLu < 3])
spp_nt_34 <- length(all_data$TLu[all_data$TLu >= 3 & all_data$TLu < 4])
spp_nt_45 <- length(all_data$TLu[all_data$TLu >= 4 & all_data$TLu < 5])
spp_nt_56 <- length(all_data$TLu[all_data$TLu < 6])

pd <- position_dodge(0.1) # move them .05 to the left and right

IS_TL <- ggplot(all_data, aes(x = reorder(TrophicSpecies, TLu), y = IS)) + 
  geom_errorbar(aes(ymin = IS-AllStrength_Q1, ymax = IS+AllStrength_Q3), width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  scale_y_log10() +
  geom_vline(xintercept = c(spp_nt_12+1, spp_nt_12+spp_nt_23+1, spp_nt_12+spp_nt_23+spp_nt_34+1, 
                            spp_nt_12+spp_nt_23+spp_nt_34+spp_nt_45+1, spp_nt_56+1), 
             linetype = "longdash", colour = "red") +
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


## Plot IS by Degree

IS_Deg <- ggplot(all_data, aes(x = reorder(TrophicSpecies, Degree), y = IS)) + 
  geom_errorbar(aes(ymin = IS-AllStrength_Q1, ymax = IS+AllStrength_Q3), width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  scale_y_log10() +
  labs(x = "Trophic species (ascending degree)", y = "Interaction strength (log scale)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15))
IS_Deg


## Interactive plot (IS, Degree, TL)

plot_3d <- plot_ly(data=all_data, x=~Degree, y=~TLu, 
                   z=~IS, type="scatter3d", mode="markers", 
                   color="black", marker = list(size = 5), hoverinfo = "text",
                   text = ~paste(TrophicSpecies, '<br>Degree:', Degree, '<br>Trophic level:', TLu, 
                                 '<br>Mean IS:', round(IS,4))) %>%  
  layout(scene = list(xaxis=list(title = "Degree"),
                      yaxis=list(title = "Trophic level"),
                      zaxis=list(title = "Interaction strength", color = "red")))
plot_3d


## Plot IS by Trophic Similarity

IS_TS <- ggplot(all_data, aes(x = reorder(TrophicSpecies, -meanTrophicSimil), y = IS)) +
  geom_errorbar(aes(ymin = IS-AllStrength_Q1, ymax = IS+AllStrength_Q3), width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  scale_y_log10() +
  labs(x = "Trophic species (descending TS)", y = "Interaction strength (log scale)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15))
IS_TS


## Plot TL by TS

TL_TS <- ggplot(all_data, aes(x = meanTrophicSimil, y = TLu)) +
  geom_point() +
  labs(x = "Trophic similarity", y = "Trophic level") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))
TL_TS

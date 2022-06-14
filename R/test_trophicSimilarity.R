require(cheddar)

cc <- LoadCommunity("Data")
cc
ts <- TrophicSimilarity(cc)
ts["Bodo saltans","Bodo saltans"]

sum(ts[,"Orcinus orca"])
vv <- ts[,"Orcinus orca"]
vv[order(desc(vv))]

mean(ts["Bodo saltans",])
mean(ts["Orcinus orca",])

#
# Genera data.frame con nombre de especie y meanTrophicSimil
#
mts <- tibble(TrophicSpecies=rownames(ts), meanTrophicSimil=colMeans(ts))
mts %>% filter(species=="Orcinus orca")

#
# Join with species traits data frame 'spp_attr'
#
spp_all <- spp_attr %>% 
  left_join(mts)
write.csv(spp_all, file = "Data/spp_all_attr.csv")

#
# Plot total interaction strength by TS
#
(plot_totalstr_ts <- ggplot(spp_all, aes(x = reorder(TrophicSpecies, -meanTrophicSimil), y = TotalStrength)) +
    geom_point() +
    labs(x = "Trophic species (descending TS)", y = "Total interaction strength") +
    # ylim(c(0,0.07)) +
    scale_y_log10() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15)))


require(plotly)

ggplotly(plot_totalstr_ts)


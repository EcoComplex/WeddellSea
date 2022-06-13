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
mts <- tibble(meanTrophicSimil=colMeans(ts),species=rownames(ts))

mts %>% filter(species=="Orcinus orca")

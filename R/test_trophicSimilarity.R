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


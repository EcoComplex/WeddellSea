#
# Functions from S. Kortsch
#

# weighted trophic level and omnvory
#
TroLev <- function(fw)
{
  fw <- t(fw)
  nn <- rowSums(fw); nn[nn==0] <- 1
  ww <- diag(1/nn)
  L1 <- ww %*% fw
  L2 <- L1 - diag(rep(1,length(nn)))
  b <- -1*rep(1,length(nn))
  Tro.lev <- solve(L2) %*% b
 
  return(Tro.lev)
}


Omnivory.species = function(i, fw, TL){
  # computes omnivory of species i
  # TL: vector of all species' TLs
  if (TL[i] == 1) {
    omn = 0
    return(omn)
    }
  prey = fw[,i] != 0
  omn = 1/sum(fw[,i]) * sum(fw[prey,i] * (TL[prey] - (TL[i] - 1))^2)
  return(omn)
}


Omnivory.web = function(fw, TL){
  # omns: omnivory value of each species in the web
  omns = sapply(1:nrow(fw), Omnivory.species, fw, TL, simplify = TRUE)
  return(omns)
}

# # fw of 5 species, sp 4 eats 2 and 3, sp 5 eats on 1 and 4
# # 1, 2, 3 are basal

# fw = matrix(0, ncol = 5, nrow = 5)
# fw[c(2,3), 4] = 1
# fw[c(2,3), 4] = 1
# fw[c(1,4), 5] = 1
# TL = TroLev(fw)
# Omnivory.web(fw, TL)
#
# # same fw as before, but 5 forage much more on 4 than on 1
# fw2 = fw
# fw2[4, 5] = 3
# TL2 = TroLev(fw2)
# Omnivory.web(fw2, TL2)

 
#function to calculate weighted indices
#weighted connectance, generality and vulnerability and the SDs of these

fluxind<- function(fluxes){
  res<-c()
# The flux matrix
W.net<-as.matrix(fluxes) #fluxmatrix from fluxweb

### Taxon-specific Shannon indices of inflows
# sum of k species inflows --> colsums
sum.in<-apply(W.net, 2, sum)

# Diversity of k species inflows
H.in.mat<- t(t(W.net)/sum.in)*t(log(t(W.net)/sum.in)) #columns divided by the total col sum
H.in.mat[!is.finite(H.in.mat)] <- 0 #converts NaN to 0's
H.in<- apply(H.in.mat, 2, sum)*-1

### Taxon-specific Shannon indices of outflows
# sum of k speies outflows --> rowsums
sum.out<-apply(W.net, 1, sum)

# Diversity of k species outflows
H.out.mat<- (W.net/sum.out)*log(W.net/sum.out) #rows divided by the total row sum
H.out.mat[!is.finite(H.out.mat)] <- 0 #converts NaN to 0's
H.out<- apply(H.out.mat, 1, sum)*-1

# Effective number of prey or resources = N(R,k) 
# The reciprocal of H(R,k) --> N (R,k) is the equivalent number of prey for species k
# N.res<-exp(H.in)
N.res<-ifelse(sum.in==0, H.in, exp(H.in))

# Effective number of predators or consumers = N(C,k) 
# The reciprocal of H(C,k) --> N (C,k) is the equivalent number of predators for species k
# N.con<-exp(H.out)
N.con<-ifelse(sum.out==0, H.out, exp(H.out))

### Quantitative Weighted and unweighted link density and weighted connectance
no.species<-ncol(W.net)
# The unweighted link density (LDw) is:
# LD.uw<- sum(N.res/no.species) + sum(N.con/no.species)/2
#I think parenthesis are missing in above formula
#I decided to follow the equation in the manuscript
res$qLD.uw<- 1/(2*no.species)* (sum(N.res) + sum(N.con))

#unweighted connectance
res$qC.uw<-res$qLD.uw/no.species

# The weighted link density (LDw) is:
# In the weighted version the effective number of predators for species i is weighted by i's 
# contribution to the total outflow
# the same is the case for the inflows

tot.mat<- sum(W.net)
# LD.w <- (sum((sum.in/tot.mat)*N.res) + sum((sum.out/tot.mat)*N.con))/2
# equivalent to next formula, but next one is closer to manuscript
res$qLD.w <- 1/(2*tot.mat)*(sum(sum.in*N.res) + sum(sum.out*N.con))

#Weighted connectance
res$qC.w<- res$qLD.w/no.species

#total number of qualitative links
res$qL <- no.species * res$qLD.uw

# positional.index
pos.ind<- sum.in*N.res/(sum.in*N.res+sum.out*N.con) #postional index
basal.sp<-pos.ind[pos.ind==0] #basal species = 0
top.sp<- pos.ind[pos.ind==1] #top predator =1
#defintion according to Bersier et al. 2002 top species = [0.99, 1]
#int.sp<- pos.ind[pos.ind>0&pos.ind<0.99]# intermediate are all species that are not basal nor top

con.sp<-length(pos.ind)-length(basal.sp)# all consumer taxa except basal
# unweighted quantitative Generality
res$qG.uw<- sum(N.res)/con.sp
# weighted quantitative Generality
res$qG.w<-sum(sum.in*N.res/sum(W.net))

res.sp<- length(pos.ind)-length(top.sp)
#unweighted quantitative Vulnerability
res$qV.uw<- sum(N.con)/res.sp
# weighted quantitative Vulnerability
res$qV.w<-sum(sum.out*N.con/sum(W.net))

#Standard deviation of stanardized unweighted quantitative Generality
s<-length(pos.ind)
res$qGsd.uw<-sd(s*N.res/sum(N.res)) #the mean  of (s*N.res/sum(N.res)  is = 1

#Standard deviation of stanardized weighted quantitative Generality
res$qGsd.w<-sd(s*sum.in*N.res/sum(N.res*sum.in))

#Standard deviation of stanardized unweighted quantitative Vulnerability
res$qVsd.uw<-sd(s*N.con/sum(N.con)) #the mean of  s*N.con/sum(N.con) is = 1

#Standard deviation of stanardized weighted quantitative Vulnerability
res$qVsd.w<-sd(s*sum.out*N.con/sum(N.con*sum.out))

return(res)
}

 
#
# END Functions from S. Kortsch
#
 

# Plot modules vs trophic level
#
#
plot_modules_TL <- function(redl,modulos)
{
  #
  # Calculate Trophic Position
  #
  require(NetIndices)
  troph.net2<-TrophInd(get.adjacency(redl,sparse=F))
  layout.matrix.1<-matrix(
    nrow=length(V(redl)),  # Rows equal to the number of vertices
    ncol=2
  )
  
  # 
  # Add colors to nodes 
  #
  require(RColorBrewer)
  
  colTL <-as.numeric(cut(troph.net2$TL,11))
  colnet <- brewer.pal(11,"RdYlGn")
  V(redl)$color <- colnet[12-colTL]
  
  #
  # Plot modules
  #
  layout.matrix.1[,2]<-jitter(troph.net2$TL,0.1) # y-axis value based on trophic level
  layout.matrix.1[,1]<-jitter(modulos$membership,1) # randomly assign along x-axis
  
#  png("Figures/RegionalFoodWeb.png",width=6,height=6,units="in",res=600)
#  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
  
  plot(redl, vertex.color=vertex_attr(redl)$cor, #vertex.label=NA,
       #vertex.size=log(3*igraph::degree(redl)),
       edge.width=.3,edge.arrow.size=.4, 
       edge.color="grey50",
       edge.curved=0.3, layout=layout.matrix.1)
  
}


#' Calculate motif counts for observed network and CI for erdos-renyi random networks and Z-scores 
#'
#' @param red igraph network object
#' @param nsim number of simulation to calculate random networks with the same nodes and links
#'
#' @return data.frame with all the results
#' @export
#'
#' @examples
calc_motif_random <- function(red, nsim=1000)
{
    Size <- vcount(red)
    Links <- ecount(red)
    
    redes.r <- lapply(1:nsim, function (x) {
      e <- erdos.renyi.game(Size, Links, type="gnm",directed = TRUE)
      while(components(e)$no>1)
        e <- erdos.renyi.game(Size, Links, type="gnm",directed = TRUE)
      
      return(e) }
    )
    
    ind <- data.frame()
    require(doParallel)
    cn <-detectCores()
    #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
    cl <- makeCluster(cn)
    registerDoParallel(cl)
    
    ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
    {
      mot <- triad_census(redes.r[[i]])
      mot[4] # Exploitative competition
      mot[5] # Apparent competition
      mot[6] # Tri-trophic chain
      mot[9] # Omnivory

     data.frame(explComp=mot[4],apprComp=mot[5],triTroph=mot[6],omnivory=mot[9])
    }
    stopCluster(cl)
    # 99% confidence interval
    #
    
    qEC <- quantile(ind$explComp,c(0.005,0.995))
    qAC <- quantile(ind$apprComp,c(0.005,0.995))
    qTT <- quantile(ind$triTroph,c(0.005,0.995))
    qOM <- quantile(ind$omnivory,c(0.005,0.995))
    
    # Calculate motif for the original network
    obs <- triad_census(red)
    
    zEC <- (obs[4] - mean(ind$explComp))/sd(ind$explComp)
    zAC <- (obs[5] - mean(ind$apprComp))/sd(ind$apprComp)
    zTT <- (obs[6] - mean(ind$triTroph))/sd(ind$triTroph)
    zOM <- (obs[9] - mean(ind$omnivory))/sd(ind$omnivory)
    
    return(data_frame(explComp=obs[4],apprComp=obs[5],triTroph=obs[6],omnivory=obs[9],zEC=zEC,zAC=zAC,zTT=zTT,zOM=zOM,EClow=qEC[1],EChigh=qEC[2],AClow=qAC[1],AChigh=qAC[2],TTlow=qTT[1],TThigh=qTT[2],OMlow=qOM[1],OMhigh=qOM[2]))         
}    

# Calculation of the clustering coefficients and average path for random network simulations
#
#
calc_modularity_random<- function(red, nsim=1000){
  
  t <- calc_topological_indices(red)
    
  redes.r <- lapply(1:nsim, function (x) {
    e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
    while(components(e)$no>1)
      e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
    
    return(e) }
    )

  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
#  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
    {
    m<-cluster_spinglass(redes.r[[i]])
    modl <- m$modularity
    ngrp <- length(m$csize)
    clus.coef <- transitivity(redes.r[[i]], type="Global")
    cha.path  <- average.path.length(redes.r[[i]])
    data.frame(modularity=modl,ngroups=ngrp,clus.coef=clus.coef,cha.path=cha.path)
  }
  stopCluster(cl)
  ind <- ind %>% mutate(gamma=t$Clustering/clus.coef,lambda=t$PathLength/cha.path,SWness=gamma/lambda)
  # 99% confidence interval
  #
  qSW <- quantile(ind$SWness,c(0.005,0.995))
  qmo <- quantile(ind$modularity,c(0.005,0.995))
  qgr <- quantile(ind$ngroups,c(0.005,0.995))
  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)
  mmo <- mean(ind$modularity)
  mgr <- mean(ind$ngroups)
  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2  
  return(list(su=data_frame(rndCC=mcc,rndCP=mcp,rndMO=mmo,rndGR=mgr,SWness=mSW,SWnessCI=mCI,MOlow=qmo[1],MOhigh=qmo[2],
                    GRlow=qgr[1],GRhigh=qgr[2]), sim=ind))         
}


#' Calc incoherence z-score and confidence interval under a random Erdos-Renyi 
#' networks with the condition of at least one basal node
#'
#' @param g igraph network object
#' @param ti trophic level vector
#'
#' @return
#' @export
#'
#' @examples
calc_incoherence_z <- function(g,ti=NULL,nsim=1000) {

    t <- calc_topological_indices(g)
    
    redes.r <- lapply(1:nsim, function (x) {
        e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
        basal <- length(V(e)[degree(e,mode="in")==0])
        while(components(e)$no>1 | basal==0){
          e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
          basal <- length(V(e)[degree(e,mode="in")==0])
        }
      return(e) }
    )
    
    ind <- data.frame()
    require(doParallel)
    cn <-detectCores()
    # #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
    cl <- makeCluster(cn)
    registerDoParallel(cl)
    
    ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph',.export = 'calc_incoherence') %do% 
    {
      m<-calc_incoherence(redes.r[[i]])
      data.frame(Q=m$Q,mTI=m$mTI)
    }
    stopCluster(cl)
    
    qQ <- quantile(ind$Q,c(0.005,0.995))
    qTI <- quantile(ind$mTI,c(0.005,0.995))
    rndQ <- mean(ind$Q)
    rndTI <- mean(ind$mTI)

    m <- calc_incoherence(g,ti)
    
    zQ <-  (m$Q- rndQ)/sd(ind$Q)
    zTI <- (m$mTI - rndTI)/sd(ind$mTI) # the same as sd(ind$mTI)
    #
    return(data_frame(rndQ=rndQ,rndTI=rndTI,Qlow=qQ[1],Qhigh=qQ[2],
                      TIlow=qTI[1],TIhigh=qTI[2],zQ=zQ,zTI=zTI))         
    
    
}
  

#' Calc topological roles 
#'
#' @param g an Igraph object with the network 
#' @param nsim  number of simulations with different community
#'
#' @return a  data frame with two fields: within_module_degree, among_module_conn
#' @export
#'
#' @examples
calc_topological_roles <- function(g,nsim=1000)
{
  
  toRol <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  toRol <- foreach(idx=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
  {
    # within-module degree
    #
    # Standarized Within module degree z-score 
    #
    m<-cluster_spinglass(g)
    spingB.mem<- m$membership
    
    l<-vector()
    memMod<-vector()
    
    for (i in 1:vcount(g)){
      
      sp.in.deg <- V(g)[nei(i, "in")]
      sp.out.deg<- V(g)[nei(i, "out")]
      mem.sp<-spingB.mem[i]
      k<- length(which(spingB.mem[c(sp.in.deg, sp.out.deg)]==mem.sp))
      mem<- which(spingB.mem==mem.sp)
      
      for (m in 1:length(mem)){
        mem.in.deg <- V(g)[nei(mem[m], "in")]
        mem.out.deg<- V(g)[nei(mem[m], "out")]
        memMod.id<- length(which(spingB.mem[c(mem.in.deg, mem.out.deg)]==mem.sp))
        memMod[m]<- memMod.id
      }
      
      k.ave<- mean(memMod)
      k.sd<- sd(memMod)
      l[i]<- (k-k.ave)/k.sd
    }
    
    # among module connectivity
    r<- vector()
    for (i in 1:vcount(g)){
      
      d<-degree(g)[i]
      sp.in.deg <- V(g)[nei(i, "in")]
      sp.out.deg<- V(g)[nei(i, "out")]
      mod.sp<-table(spingB.mem[c(sp.in.deg, sp.out.deg)])
      mod.no<-as.numeric(names(mod.sp))
      mod<-rep(0, length(unique(spingB.mem)))
      mod[mod.no]<-c(mod.sp)
      r[i]<- 1-(sum((mod/d)^2))
    }
    return(data.frame(node=1:vcount(g),within_module_degree=l, among_module_conn=r))
    
  }
  stopCluster(cl)
  
  # toRol %>% group_by(node) %>% summarize(wtmLowCI=quantile(within_module_degree,0.005,na.rm=TRUE),
  #                                        wtmHiCI=quantile(within_module_degree,0.995,na.rm=TRUE),
  #                                        amcLowCI=quantile(among_module_conn,0.005,na.rm=TRUE),
  #                                        amcHiCI=quantile(among_module_conn,0.995,na.rm=TRUE),
  #                                        within_module_degree=mean(within_module_degree,na.rm=TRUE),
  #                                        among_module_conn=mean(among_module_conn,na.rm=TRUE))
  return(toRol)
}


#' Plot topological roles
#'
#' @param tRoles Calculated topological roles with the function calc_topological_roles
#' @param g Igraph network object   
#' @param spingB Igraph community object
#'
#' @return
#' @export
#'
#' @examples
plot_topological_roles <- function(tRoles,g,spingB){
  
  
  spingB.mem<- spingB$membership
  
  l <- tRoles$within_module_degree
  r <- tRoles$among_module_conn
  # Plot
  require(RColorBrewer)
  colbar.FG <- brewer.pal(length(spingB$csize),"Dark2")
  groups.B<- spingB.mem # 1=phytopl, 2=zoopl, 3=benthos, 4=fish, 5=birds, 6=mammals
  
  par(mfrow=c(1,1))
  par(oma=c(2.7,1.3,0.7,1)) # change outer figure margins
  par(mai=c(1,1,0.7,0.7)) # size of figure margins
  yran <- range(l)
  xran <- range(r)
  plot(l~r, type="n", axes=T ,tck=F, lwd=2, ann=F, cex.axis=1.2, xlim=xran, ylim=yran)
  lines(c(0.625,0.625), yran, col="grey")
  lines(xran, c(2.5, 2.5), col="grey")
  points(r, l, col=colbar.FG[groups.B], pch=20, cex=2)
  mtext(2, text="Within module degree", line=4,font=2, cex=1.2)
  mtext(1, text="Among module connectivity",line=4, font=2, cex=1.2)
  axis(1, tck=0.02, lwd=1,las=2,lty=1, labels=F, xlim=c(0,1))
  axis(2,tck=0.02, labels=FALSE)
  
  # Which are the module hubs: many links within its own module.
  #
  modhub <- which(l>2.5)
  modhub <- modhub[which(l>2.5) %in% which(r<=0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  hub_conn <- data.frame()
  
  if(length(modhub)) {
    text(r[modhub],l[modhub],labels = modlbl,cex=0.7,pos=3)
    hub_conn <- data.frame(type="modhub",node=modhub,name=modlbl)  
  }
  
  #points(r[modhub], l[modhub], cex=4, col="blue", pch=20)
  
  # Which are the hub connectors: high within and between-module connectivity
  #                              and are classified super-generalists
  #
  modhub <- which(l>2.5)
  modhub <- modhub[which(l>2.5) %in% which(r>0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  if(length(modhub)) {
    text(r[modhub],l[modhub],labels = modlbl,cex=0.7,pos=3)
  }
  
  #points(r[modhub], l[modhub], cex=4, col="blue", pch=20)
  if(length(modhub)){
    hub_conn <- rbind(hub_conn, data.frame(type="hubcon",node=modhub,name=modlbl))  
  }
  
  # Which are the module specialist: Few links and most of them within its own module
  #
  modhub <- which(l<=2.5)
  modhub <- modhub[which(l<=2.5) %in% which(r<=0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  
  hub_conn <- rbind(hub_conn, data.frame(type="modspe",node=modhub,name=modlbl))  
  
  # Which are the module connectors: Few links and between modules
  #
  modhub <- which(l<=2.5)
  modhub <- modhub[which(l<=2.5) %in% which(r>0.625)]
  modlbl <- unlist(vertex_attr(g,index=modhub))
  if(is.null(modlbl))
    modlbl <- modhub
  
  hub_conn <- rbind(hub_conn, data.frame(type="modcon",node=modhub,name=modlbl))  
  
}


#' Calculate average topological roles doing nsimStep simulations and repeating until there is no
#' differences using an Anderson Darling test.
#'
#' @param g igraph object with the network 
#' @param nsimStep number of repeated simulations until testing for differences, the minimun number of simulations is nsimStep*2
#'
#' @return data.frame with topological roles averaged over n*nsimStep repetitions
#' @export
#'
#' @examples
calc_avg_topological_roles <- function(g, net_name,nsimStep){
  
  tR1 <- calc_topological_roles(g,nsimStep)                       # 30 simulations is enough to obtain stable topological roles
  tsim <- nsimStep
  topoRoles_mWA_temp  <- tR1 %>% group_by(node) %>% summarize(wtmLowCI=quantile(within_module_degree,0.005,na.rm=TRUE),
                                              wtmHiCI=quantile(within_module_degree,0.995,na.rm=TRUE),
                                              amcLowCI=quantile(among_module_conn,0.005,na.rm=TRUE),
                                              amcHiCI=quantile(among_module_conn,0.995,na.rm=TRUE),
                                              within_module_degree=mean(within_module_degree,na.rm=TRUE),
                                              among_module_conn=mean(among_module_conn,na.rm=TRUE))
  
  
  print(tsim)
  # Loop
  while(TRUE){
    tR1 <- bind_rows(tR1, calc_topological_roles(g,nsimStep))
    
    saveRDS(tR1,"TopoRoles_mWA_temp.rds")
    
    #tR1 <- readRDS("TopoRoles_mWA_temp.rds")
    
    tR  <- tR1 %>% group_by(node) %>% summarize(wtmLowCI=quantile(within_module_degree,0.005,na.rm=TRUE),
                                                wtmHiCI=quantile(within_module_degree,0.995,na.rm=TRUE),
                                                amcLowCI=quantile(among_module_conn,0.005,na.rm=TRUE),
                                                amcHiCI=quantile(among_module_conn,0.995,na.rm=TRUE),
                                                within_module_degree=mean(within_module_degree,na.rm=TRUE),
                                                among_module_conn=mean(among_module_conn,na.rm=TRUE))
    
    
    require(kSamples)
    t1 <- ad.test(list(tR$among_module_conn,topoRoles_mWA_temp$among_module_conn),method="simulated",nsim=1000)
    t2 <- ad.test(list(tR$within_module_degree,topoRoles_mWA_temp$within_module_degree),method="simulated",nsim=1000)
    topoRoles_mWA_temp <- tR %>% mutate(Network=net_name)
    tsim <- tsim +nsimStep
    print(tsim)
    if(t1$ad[1,4]>0.1 && t2$ad[1,4]>0.1) break()
    
  }
  
  return(topoRoles_mWA_temp)
}

#' Plot topological roles by trophic level and module
#'
#' @param netFrame dataframe with all the networks 
#' @param netName String with name of the food web to analyse
#' @param deadNodes Vector of strings with name of dead nodes to calculate trophic level
#' @param modulObj Igraph community object with the module organization of the food web
#' @param topoFrame dataframe with topological role and node index
#' @param legendPos position of the legend "topleft", "topright" or if "" no legend.
#' @param redl igraph object, if it is not null the network is taken from it    
#'
#' @return
#' @export
#'
#' @examples

plotTopoRolesByTLByMod <- function(netFrame,netName,deadNodes,modulObj,topoFrame,legendPos="",redl=NULL){
  # 
  # Local 
  #
  if(is.null(redl)){
    
    dtot1 <- as.matrix(netFrame %>% filter(Network==netName) %>% dplyr::select(Prey_name,Predator_name))
    redl <- graph_from_edgelist(dtot1, directed  = T)
    redl <- simplify(redl)
  }
  require(NetIndices)
  if(deadNodes=="")
      troph.net2<-TrophInd(get.adjacency(redl,sparse=F))
  else
      troph.net2<-TrophInd(get.adjacency(redl,sparse=F),Dead=deadNodes)
  layout.matrix.1<-matrix(
    nrow=length(V(redl)),  # Rows equal to the number of vertices
    ncol=2
  )
  
  # 
  # Add colors with topological roles to nodes 
  #
  require(RColorBrewer)
  colnet <- brewer.pal(4,"Paired")
  
  hc <- topoFrame %>% mutate(type = factor(type)) %>% filter(Network==netName) %>% arrange(node) %>% mutate(col= as.numeric(type), TL=troph.net2[,1]) 
  V(redl)$color <- colnet[hc$col]
  
  # Transform y-axis coordinates
  #
  maxnew <- max(hc$TL)
  minnew <- min(hc$TL)
  maxold <- 1
  minold <- -1
  t2 <- function(x) (maxold-minold)/(maxnew -minnew)*(x - maxnew)+maxold 
  
  
  #
  # Plot modules
  #
  layout.matrix.1[,2]<-jitter(troph.net2$TL,0.4) # y-axis value based on trophic level
  layout.matrix.1[,1]<-jitter(modulObj$membership,1) # randomly assign along x-axis
  
  
  plot(redl, vertex.color=vertex_attr(redl)$cor,vertex.label=NA,
       vertex.size=log(3*igraph::degree(redl)),
       edge.width=.3,edge.arrow.size=.2, 
       edge.color=add.alpha("grey50",0.5),
       edge.curved=0.3, layout=layout.matrix.1)
  
  
  axis(side=2,at=t2(1:4),labels=1:4,  las=1, col = NA, col.ticks = 1)
  
  legstr <- levels(hc$type)
  legstr <- c("Hub conn.", "Mod. Conn.", "Mod. Hubs", "Mod. Spec.")
  if(legendPos!="")
      legend(legendPos, pch=19, col=colnet, legend= legstr)
  
}


# Add alpha to base plot colors
# 
#
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}




getTopoRolesTLdegree <- function(netFrame,netName,deadNodes,topoFrame,topoType=NULL){
  # 
  # Igraph object from dataframe
  #
  dtot1 <- as.matrix(netFrame %>% filter(Network==netName) %>% dplyr::select(Prey_name,Predator_name))
  redl <- graph_from_edgelist(dtot1, directed  = T)
  redl <- simplify(redl)
  
  require(NetIndices)
  TL<-TrophInd(get.adjacency(redl,sparse=F),Dead=deadNodes)

  if(!is.null(topoType)){
    topoFrame %>% filter(type==topoType,Network==netName) %>% rowwise() %>% mutate( preys=degree(redl,node,mode=c("in")), predators= degree(redl,node,mode=c("out")), trophLevel=TL[node,1])
  } else {
    topoFrame %>% filter(Network==netName) %>% group_by(type) %>% rowwise() %>% mutate( preys=degree(redl,node,mode=c("in")), predators= degree(redl,node,mode=c("out")), trophLevel=TL[node,1])
  }
}
  
#' Title Plot net assembly model time series of S and L only the last timeW steps are ploted
#'
#' @param AA output of a net assembly model
#' @param timeW time window used
#' @param fname file name to save the plot 
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel <- function(AA,timeW,fname=NULL){
  require(viridis)
  colnet <- viridis(3)
  
  tf <- length(AA$L)
  if(tf<timeW) stop("timeW parameter must be less than the time of the simulation")
  
  dfA <- data.frame(S=AA$S[(tf-timeW):tf],L=as.numeric(AA$L[(tf-timeW):tf]),T=c((tf-timeW):tf))
  dfA$C <- dfA$L/(dfA$S*dfA$S)
  if(is.null(fname)){
    gS <- ggplot(dfA, aes(x=T,y=S)) + geom_line(colour=colnet[1]) + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype = 2,colour="grey50")
    print(gS)
    gL <- ggplot(dfA, aes(x=T,y=L)) + geom_line(colour=colnet[2]) + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype = 2,colour="grey50")
    print(gL)
    gC <- ggplot(dfA, aes(x=T,y=C)) + geom_line(colour=colnet[3]) + theme_bw() + ylab("C") + geom_hline(yintercept = mean(dfA$C),linetype = 2,colour="grey50")
    print(gC)
    return(list(gS=gS,gL=gL,gC=gC))
  } else {
    require(cowplot)
    g1 <- ggplot(dfA, aes(x=T,y=S)) + geom_line() + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype=3)
    g2 <- ggplot(dfA, aes(x=T,y=L)) + geom_line() + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype=3)
    g3 <- plot_grid(g1,g2,labels = c("A","B"),align = "h")
    save_plot(fname,g3,base_width=8,base_height=5,dpi=600)
  }
}


#' Title Plot net assembly model S and L average by a moving window to check if equilibrium is reached
#'
#' @param AA output of a net assembly model
#' @param timeW time window used
#' @param fname file name to save the plot 
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel_eqw <- function(AA,timeW,fname=NULL){

  df <- data.frame(S=AA$S,L=as.numeric(AA$L),T=c(1:tf))
  grandS <- mean(df$S[timeW:nrow(df)])
  grandL <- mean(df$L[timeW:nrow(df)])
  
  df$gr <- rep(1:(nrow(df)/timeW), each = timeW)
  df <- df %>% group_by(gr) %>% summarise(mS=mean(S),sdS=sd(S), mL=mean(L), sdL=sd(L),time=max(T))
  if(is.null(fname)){
    print(ggplot(df,aes(y=mS,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mS-sdS,ymax=mS+sdS)) + scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandS,linetype=3 ))
    print(ggplot(df,aes(y=mL,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mL-sdL,ymax=mL+sdL))+ scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandL,linetype=3 ))
  } else {
    require(cowplot)
    g1 <- ggplot(df,aes(y=mS,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mS-sdS,ymax=mS+sdS)) + scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandS,linetype=3 )
    g2 <- ggplot(df,aes(y=mL,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mL-sdL,ymax=mL+sdL))+ scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandL,linetype=3 )
    g3 <- plot_grid(g1,g2,labels = c("A","B"),align = "h")
    save_plot(fname,g3,base_width=8,base_height=5,dpi=600)
  }
    
  return(df)
}



#' Estimation of z-scores using Meta-Web assembly model as a null 
#'
#' @param red This is the reference network as an igraph object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param ti  trophic level vector 
#' @param nsim number of simulations
#'
#' @return
#' @export
#'
#' @examples
calc_modularity_metaWebAssembly<- function(red, Adj, mig,ext,nsim=1000,ti=NULL){
  
  t <- calc_topological_indices(red)
  final_time <- 500  # Final time used in simulations of the meta-web assembly
  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))
  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages=c('MetaWebAssemblyModels','igraph'), 
                 .export = c('Adj','ext','mig','final_time','calc_incoherence')) %dopar% 
  {
    AA <- metaWebNetAssembly(Adj,mig,ext,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    # Select only a connected subgraph graph 
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    mmm<-cluster_spinglass(g)
    modl <- mmm$modularity
    ngrp <- length(mmm$csize)
    clus.coef <- transitivity(g, type="Global")
    cha.path  <- average.path.length(g)
    mmm<-calc_incoherence(g)

    data.frame(modularity=modl,ngroups=ngrp,clus.coef=clus.coef,cha.path=cha.path,Q=mmm$Q,mTI=mmm$mTI)
  }
  stopCluster(cl)
  ind <- ind %>% mutate(gamma=t$Clustering/clus.coef,lambda=t$PathLength/cha.path,SWness=gamma/lambda)
  # 99% confidence interval
  #
  qSW <- quantile(ind$SWness,c(0.005,0.995))
  qmo <- quantile(ind$modularity,c(0.005,0.995))
  qgr <- quantile(ind$ngroups,c(0.005,0.995))
  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)
  mmo <- mean(ind$modularity)
  mgr <- mean(ind$ngroups)
  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2  

  qQ <- quantile(ind$Q,c(0.005,0.995))
  qTI <- quantile(ind$mTI,c(0.005,0.995))
  mdlQ <- mean(ind$Q)
  mdlTI <- mean(ind$mTI)
  
  m <- calc_incoherence(red,ti)
  
  zQ <-  (m$Q- mdlQ)/sd(ind$Q)
  zTI <- (m$mTI - mdlTI)/sd(ind$mTI) # the same as sd(ind$mTI)
  
  return(list(su=data_frame(mdlCC=mcc,mdlCP=mcp,mdlMO=mmo,mdlGR=mgr,SWness=mSW,SWnessCI=mCI,MOlow=qmo[1],MOhigh=qmo[2],
                    GRlow=qgr[1],GRhigh=qgr[2], mdlQ=mdlQ,mdlTI=mdlTI,Qlow=qQ[1],Qhigh=qQ[2],
                                                 TIlow=qTI[1],TIhigh=qTI[2],zQ=zQ,zTI=zTI,MOsd=sd(ind$modularity)),sim=ind))         
}



#' Estimation of QSS z-scores using Meta-Web assembly model as a null 
#'
#' @param red This is the reference network as an igraph object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param nsim number of simulations
#'
#' @return
#' @export
#'
#' @examples
calc_qss_metaWebAssembly<- function(red, Adj, mig,ext,nsim=1000,ncores=0){
  
  t <- calc_QSS(red,10000,ncores)
  final_time <- 500  # Final time used in simulations of the meta-web assembly
  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))
  ind <- data.frame()

  require(doFuture)
  registerDoFuture()
  if(ncores) {
    cn <- future::availableCores()
    if(ncores>cn)
      ncores <- cn
    plan(multiprocess, workers=ncores)
    # cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug
    # cl <- parallel::makeCluster(ncores)
    # doParallel::registerDoParallel(cl)
    # on.exit(parallel::stopCluster(cl))
  } else {
    plan(sequential)
  }

  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE) %do% 
    {
      AA <- metaWebNetAssembly(Adj,mig,ext,final_time)
      g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
      dg <- components(g)
      g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
      
      size <- vcount(g)
      links <- ecount(g)
      
      # Select only a connected subgraph graph 
      print(paste("Sim:",i, "Size:", size))
      bind_cols(data.frame(Size=size,Links=links),calc_QSS(g,10000,ncores))
    }
  # 99% confidence interval
  #
  q_qss <- quantile(ind$QSS,c(0.005,0.995),na.rm = TRUE)
  m_qss <- mean(ind$QSS)
  q_meing <- quantile(ind$MEing,c(0.005,0.995),na.rm = TRUE)
  m_meing <- mean(ind$MEing)
  
  zQSS <- (t$QSS - m_qss)/sd(ind$QSS) # the same as sd(ind$mTI)
  zMEing <- (t$MEing - m_meing)/sd(ind$MEing)
  return(list(su=data_frame(QSS=t$QSS,mdlQSS=m_qss,QSSlow=q_qss[1],QSShigh=q_qss[2],
                zQSS=zQSS,MEing=t$MEing,mdlMEing=m_meing,MEingLow=q_meing[1],MEingHigh=q_meing[2],zMEing=zMEing),sim=ind))         
}


#' Calculate motif counts for observed network and CI for meta-web assembly model networks and Z-scores 
#'
#' @param red igraph network object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param nsim number of simulation to calculate random networks with the same nodes and links
#'
#' @return data.frame with all the results
#' @export
#'
#' @examples
calc_motif_metaWebAssembly<- function(red, Adj, mig, ext, nsim=1000)
{
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)
  
  final_time <- 500  # Final time used in simulations of the meta-web assembly
  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))

  ind <- data.frame()
  require(doParallel)
  cn <-detectCores()
  #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
  cl <- makeCluster(cn)
  registerDoParallel(cl)

  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages=c('MetaWebAssemblyModels','igraph'), 
                 .export = c('Adj','ext','mig','final_time')) %dopar% 
  {
    AA <- metaWebNetAssembly(Adj,mig,ext,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    # Select only a connected subgraph graph 
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    
    mot <- triad_census(g)
    mot[4] # Exploitative competition
    mot[5] # Apparent competition
    mot[6] # Tri-trophic chain
    mot[9] # Omnivory
    
    data.frame(explComp=mot[4],apprComp=mot[5],triTroph=mot[6],omnivory=mot[9])
  }
  stopCluster(cl)
  # 99% confidence interval
  #
  
  qEC <- quantile(ind$explComp,c(0.005,0.995))
  qAC <- quantile(ind$apprComp,c(0.005,0.995))
  qTT <- quantile(ind$triTroph,c(0.005,0.995))
  qOM <- quantile(ind$omnivory,c(0.005,0.995))
  
  # Calculate motif for the original network
  obs <- triad_census(red)
  
  zEC <- (obs[4] - mean(ind$explComp))/sd(ind$explComp)
  zAC <- (obs[5] - mean(ind$apprComp))/sd(ind$apprComp)
  zTT <- (obs[6] - mean(ind$triTroph))/sd(ind$triTroph)
  zOM <- (obs[9] - mean(ind$omnivory))/sd(ind$omnivory)
  
  return(data_frame(explComp=obs[4],apprComp=obs[5],triTroph=obs[6],omnivory=obs[9],zEC=zEC,zAC=zAC,zTT=zTT,zOM=zOM,EClow=qEC[1],EChigh=qEC[2],AClow=qAC[1],AChigh=qAC[2],TTlow=qTT[1],TThigh=qTT[2],OMlow=qOM[1],OMhigh=qOM[2]))         
}    



#' Title Plot 5 simulations of net assembly model time series of S and L only the last timeW steps are ploted
#'
#' @param metaW meta-web adjacency matrix 
#' @param m     migration
#' @param q     probability of link
#' @param a     extinction
#' @param timeW time window used
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel_sims <- function(metaW,m, q, a, tf,timeW){
  require(viridis)

  if(tf<timeW) stop("timeW parameter must be less than the time of the simulation")
  
  dfA <- data.frame()
  
  for(n in 1:5){
    AA <- metaWebNetAssembly(metaW,m,q,a,tf)
    tdfA <- data.frame(S=AA$S[(tf-timeW):tf],L=as.numeric(AA$L[(tf-timeW):tf]),T=c((tf-timeW):tf))
    tdfA$C <- tdfA$L/(tdfA$S*tdfA$S)
    tdfA$sim <- n
    dfA <- bind_rows(dfA,tdfA)
  }
  gS <- ggplot(dfA, aes(x=T,y=S,colour=sim)) + geom_point() + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gS)
  gL <- ggplot(dfA, aes(x=T,y=L,colour=sim)) + geom_point() + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gL)
  gC <- ggplot(dfA, aes(x=T,y=C,colour=sim)) + geom_point() + theme_bw() + ylab("C") + geom_hline(yintercept = mean(dfA$C),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gC)
  return(list(gS=gS,gL=gL,gC=gC))
}

#' Show the names of basal species 
#'
#' @param g igraph network
#'
#' @return
#' @export
#'
#' @examples
basal_species <- function(g){
  deg <- degree(g, mode="in") # calculate the in-degree: the number of preys
  
  V(g)$indegree <-  deg
  
  basal <- V(g)[indegree==0]
  
  V(g)[basal]$name
  
}


#' Substract duplicated from wedd to meta if it do not exist in pott
#'
#' @param meta 
#' @param wedd 
#' @param pott 
#'
#' @return igraph object
#' @export
#'
#' @examples
delete_equal_edges <- function(meta,wedd,pott)
{
  s <- vcount(wedd)
  del <- 0
  meta1 <- meta
  for(i in seq_len(s)){
    nam <- V(wedd)[i]$name
    if( nam %in% V(meta1)$name ){
      li <- incident(meta1, nam)
      if( all( as_ids(incident(wedd, nam)) %in% as_ids(li) ) ) {
        if( !(nam %in% V(pott)$name ))  {
          
          cat( i, nam, "\n") 
          cat( "S:", vcount(meta1), "L:", ecount(meta1), "\n") 
          #x <- toupper(scan(what=character(),nmax=1,quiet=TRUE))
          if( !(nam %in% c("Sediment","Detritus"))){
            meta1 <- meta1 - li
            meta1 <- meta1 - nam
            cat( "DELETED S:", vcount(meta1), "L:", ecount(meta1), "\n") 
            del <- del + 1
          }
        }
      }
    } else {
      cat( "DUPLICATED:", nam , "\n") 
    }
  }
  cat("Deleted: ", del, "\n")
  return(meta1)
  
}

#' Select the parameter set by tolerance of Species and links
#'
#' @param webs 
#' @param meta 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
sel_tolerance_local_web <- function( webs, meta, tol,simMeta){
  require(dplyr)
  require(ggplot2)
  df <- lapply(webs, function(w){
    size <- vcount(w)
    links <- ecount(w)
    sim <- simMeta %>% group_by(m,a) %>% summarize(S=mean(S),L=mean(L),C=mean(C),cost = sqrt((S - size)^2)/size+sqrt((L-links)^2)/links) %>% mutate(alpha=m/a) %>% filter(L>links*(1-tol),L<links*(1+tol),S>size*(1-tol),S<size*(1+tol)) %>% arrange(cost)
    
    print(ggplot(sim,aes(S, L))+ geom_point(alpha=0.5,shape = 21,fill="white") +  geom_point(data=data.frame(size,links),aes(size, links),shape = 21, size=2) + theme_bw() + scale_color_viridis_d())
    
    res <- sim %>% ungroup() %>% summarise_if(is.numeric, mean)
    res <- bind_cols(res, sim %>% ungroup() %>% summarise(n=n()))
  })
  #
  # Selectin simulations for Potter Cove with 30% tolerance for both Species and Links
  #
  
  do.call(rbind,df)
  
}


#' Make simulations of the metaWebAssembly model 
#' 
#' Results 
#'
#' @param mini 
#' @param mend 
#' @param mstep 
#' @param aini 
#' @param aend 
#' @param astep 
#' @param nrep number of repetitions for each parameter combinations
#' @param tf   final time for simulations
#' @param A    metaweb
#'
#' @return
#' @export
#'
#' @examples
simul_metaWebAssembly <- function(mini,mend, mstep, aini, aend, astep, nrep, tf,A){
  require(dplyr)
  # Combine parameters
  #
  arguments <- expand.grid(m = seq(from=mini,to=mend,by=mstep), a = seq(from=aini,to=aend,by=astep))

  # repeat nrep times
  #
  arguments <- arguments %>% slice(rep(row_number(), nrep))
  
  sim <- data.frame()
  require(doFuture)
  registerDoFuture()
  cn <- future::availableCores()
  plan(multiprocess, workers=cn)
  
  nsim <- dim(arguments)[1]
  sim <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE) %dopar% 
    {
      mm <- rep( arguments[i,1], nrow(A))
      aa <- rep( arguments[i,2], nrow(A))
      AA <- metaWebNetAssembly(A,mm,aa,tf)
      dfA <- data.frame(S=AA$S[200:tf],L=as.numeric(AA$L[200:tf]),T=c(200:tf))
      data.frame(m=mean(mm),a=mean(aa), S=mean(dfA$S),L=mean(dfA$L),C=mean(dfA$L)/(mean(dfA$S)*mean(dfA$S)))
    }
  return(sim)
}



#' Plot networks by trophic level with some nodes highlighted 
#'
#' @param redl igraph object, if it is not null the network is taken from it    
#' @param modulObj Igraph community object with the module organization of the food web
#' @param topoFrame dataframe with topological role and node index
#' @param legendPos position of the legend "topleft", "topright" or if "" no legend.
#'
#' @return
#' @export
#'
#' @examples

plot_troph_level_mod <- function(redl,modulObj,groups,legendPos=""){
  # 
  # Local 
  #
  if(class(redl)!="igraph"){
    stop("parameter redl must be an igraph object")
  }

  require(NetIndices)
  troph.net2<-TrophInd(get.adjacency(redl,sparse=F))

  layout.matrix.1<-matrix(
    nrow=length(V(redl)),  # Rows equal to the number of vertices
    ncol=2
  )
  
  # 
  # Add colors with topological roles to nodes 
  #
  require(RColorBrewer)
  no_groups <- length(unique(groups))
  colnet <- brewer.pal(no_groups,"Paired")
  
  V(redl)$color <- colnet[groups]
  
  # Transform y-axis coordinates
  #
  maxnew <- max(hc$TL)
  minnew <- min(hc$TL)
  maxold <- 1
  minold <- -1
  t2 <- function(x) (maxold-minold)/(maxnew -minnew)*(x - maxnew)+maxold 
  
  
  #
  # Plot modules
  #
  layout.matrix.1[,2]<-jitter(troph.net2$TL,0.4) # y-axis value based on trophic level
  layout.matrix.1[,1]<-jitter(modulObj$membership,1) # randomly assign along x-axis
  
  
  plot(redl, vertex.color=vertex_attr(redl)$cor,vertex.label=NA,
       vertex.size=log(3*igraph::degree(redl)),
       edge.width=.3,edge.arrow.size=.2, 
       edge.color=add.alpha("grey50",0.5),
       edge.curved=0.3, layout=layout.matrix.1)
  
  
  axis(side=2,at=t2(1:4),labels=1:4,  las=1, col = NA, col.ticks = 1)
  
  legstr <- unique(groups)
  legstr <- c("Hub conn.", "Mod. Conn.", "Mod. Hubs", "Mod. Spec.")
  if(legendPos!="")
    legend(legendPos, pch=19, col=colnet, legend= legstr)
  
}

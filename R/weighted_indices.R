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

#weighted trophic level and omnvory
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
fw = matrix(0, ncol = 5, nrow = 5)
fw[c(2,3), 4] = 1
fw[c(2,3), 4] = 1
fw[c(1,4), 5] = 1
TL = TroLev(fw)
Omnivory.web(fw, TL)
# 
# # same fw as before, but 5 forage much more on 4 than on 1
fw2 = fw
fw2[4, 5] = 3
TL2 = TroLev(fw2)
Omnivory.web(fw2, TL2)

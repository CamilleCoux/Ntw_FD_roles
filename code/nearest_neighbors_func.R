# edit to try out calculating nearest neighbours


# Pollinators
# so here the coordinates are calculated for the whole community, it's the centroid
# that's going to change coordinates. And then I'll just need to adapt the nearest 
# neighbour part such that only the distances to *present* neighbours are calculated.
species.coords <- function(traitmat, abund, weights, ...){
  # calculating species' coordinates and the centroid
  coords <- dbFD(traitmat, abund, w=weights, corr="cailliez", print.pco=T)$x.axes
  # now first calculating the centroid of the pollis present for each site
  centr <- list(NULL)
  # centr = averaged trait values for pollis present. 14 dimensions for 21 sites.
  for(i in 1:length(abund[,1])){
    pres <- which(abund[i,]>0)
    vec <- coords[as.logical(abund[i,]),]
    w <- abund[i,pres]
    centr[[i]] <- apply(vec, 2, weighted.mean, w = w)
  }
  centr <- do.call(rbind, centr)
  rownames(centr) <- rownames(abund)
  return(list(coords=coords, centr=centr))
}


# now to calculate originality, redundancy and uniqueness

FD_measures <- function(coords, centr, abund){
  # adds a line for the centroid coodinates in each community
  coords2  <- list(NULL)
  for (i in 1:length(abund[,1])){
    coords2[[i]] <- data.frame(rbind(coords, centr[i,]))
    rownames(coords2[[i]])[dim(coords)[1]+1] <- "centr"
  }
  
  dist.poll    <- lapply(coords2, function(x){
    d.poll<- as.matrix(dist(x, diag=TRUE, upper=TRUE))
    for (i in 1:dim(d.poll)[1]) {d.poll[i,i] <- NA}
    return(d.poll)
  })
  
  # 6. Originality: Distance to centroid of the species present
  
  orig.poll<- lapply(1:length(abund[,1]), function(i){
    x <- dist.poll[[i]]
    pres <- which(abund[i,]>0)
    x[pres,dim(x)[1]]
  })
  
  # 7. Redundancy: nearest neighbour among present species
  
  redund.poll <- lapply(1:length(abund[,1]), function(i){
    pres <- which(abund[i,]>0)
    x <- orig.poll[[i]]
    y <- dist.poll[[i]][pres,pres]
    measures <- data.frame(cbind(x, redund.poll=rep(NA,length(pres))))
    colnames(measures)[1] <- "orig.poll"
    for (j in 1:length(pres)){measures[j,2] <- min(y[j,], na.rm=T)}
    return (measures)
  })
  
  # 8. Uniqueness: mean distance to the rest of the species present
  
  uniq.poll <- lapply(1:length(abund[,1]),function(i){
    pres <- which(abund[i,]>0)
    uniq <- colMeans(dist.poll[[i]][pres,pres], na.rm=T)
    return(uniq)
  })
  
  # all together
  measures <- lapply(1:length(abund[,1]), function(i){
    m <- cbind(redund.poll[[i]],uniq.poll[[i]])
    colnames(m)[3] <- "uniq.poll"
    #m <- m[-17,]
    return(m)
  })
  measures <- do.call(rbind, measures)
  return(data.frame(
    orig = measures$orig.poll, 
    redund = measures$redund.poll, 
    uniq = measures$uniq.poll))
}
# pol.coords <- species.coords(t3, a, weight)
# pol.measures <- FD_measures(pol.coords$coords, pol.coords$centr, a)
# cor(pol.measures)
# # POLLIS:
# # originality and uniqueness: 0.75
# # redundancy and uniqueness: 0.78
# # originality and redundancy: 0.30
# pl.coords <- species.coords(pl_t,p,T,pl.weight)
# pl.measures <- FD_measures(pl.coords$centr, pl.coords$centr, p)
# cor(pl.measures) 
# # PLANTS:
# # originality and uniqueness: 0.65
# # redundancy and uniqueness: 0.92
# # originality and redundancy: 0.75
# 
# # with weighted centroid: only originality should change
# w.pol.coords <- species.coords(t3,a2,weight)
# w.pol.measures <- FD_measures(w.pol.coords$coords, w.pol.coords$centr, a2)
# cor(w.pol.measures)
# # originality and redundancy: -0.22
# # originality and uniqueness: -0.05
# 
# # # this will have to wait for romina's reply
# # w.pl.coords <- species.coords(pl_t,p,pl.weight)
# # w.pl.measures <- FD_measures(pl.coords$centr, pl.coords$centr, p)
# # cor(pl.measures) 
# # 


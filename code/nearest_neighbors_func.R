#################################################################################
# code by Camille Coux
# last modified 20 Jan 2016
#
# 2 functions to calculate species' functional originality and uniqueness:
# - species.coords calls the dbFD function from FD package to calculate species'
#   coordinates in the trait space, and rearranges the outputs
# - FD_measures calculates originality and uniqueness for species present in 
#   a given community.
# 
# Inputs are the same as those needed for dbFD; row/column name-matching applies.
# Output: a dataframe with species as rows and FD measures as columns, in addition
# to abundances and farm (site name).
# 
################################################################################
################################################################################



# The coordinates are calculated for the whole community, it's the centroid
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


# now to calculate originality and uniqueness

FD_measures <- function(coords, centr, abund){
  # adds a line for the centroid coodinates in each community
  coords2  <- list(NULL)
  for (i in 1:length(abund[,1])){
    coords2[[i]] <- data.frame(rbind(coords, centr[i,]))
    rownames(coords2[[i]])[dim(coords)[1]+1] <- "centr"
  }
  
  dists_centr   <- lapply(coords2, function(x){
    d.poll<- as.matrix(dist(x, diag=TRUE, upper=TRUE))
    for (i in 1:dim(d.poll)[1]) {d.poll[i,i] <- NA}
    return(d.poll)
  })
  
  # Originality: Distance to centroid of the species present
  
  originality<- unlist(lapply(1:length(abund[,1]), function(i){
    dists <- dists_centr[[i]]
    pres <- which(abund[i,]>0)
    dists[pres,dim(dists)[1]]
  }))
  
  # Uniqueness: nearest neighbour among present species
  
  uniqueness <- unlist(lapply(1:length(abund[,1]), function(i){
    pres <- which(abund[i,]>0)
    dists <- dists_centr[[i]][pres,pres]
    uniq <- NULL
    for (j in 1:length(pres)){uniq <- c(uniq, min(dists[j,], na.rm=T))}
    return (uniq)
  }))
  
  # also need to keep track of farm, species names and their abundances
  
  measures <- do.call(rbind, lapply(1:length(abund[,1]), function(i){
    pres <- which(abund[i,]>0)
    abundance <- abund[i, pres]
    sp <- names(abundance)
    site <- rownames(abund)[i]
    farm <-  rep(site, length(pres))
    tab <- data.frame(farm=farm, sp=sp, abundance=as.numeric(abundance))
    return (tab)
  }))
  
  # merge all
  measures$orig <- originality
  measures$uniq <- uniqueness
  
  return(as.data.frame(measures))
  
}

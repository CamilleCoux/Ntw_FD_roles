# code by Camille Coux
# last modified 21 Dec 2015

# FD and bipartite metric calculation with ntw weighed by the number of insect 
# individuals from a same species carrying pollen of 1 plant. 


library(FD)
library(ESM)
library(bipartite)
library(magrittr)

source("nearest_neighbors_func.R")

# 1. Import data:
# poll traits
t3 <- read.table("pollinator_traits.csv", sep=",", header=T, row.names = 1 )

# assigning polli trait weights
weight <- c(0.5, 0.5,1,0.5,0.5,1, rep(1/6, 6), rep(1,3))

# plant traits
pl_t <-  read.table("plant_traits.csv", sep=",", header=T, row.names = 1 )

# assigning plant trait weights
pl.weight <- c(rep(1, 8), rep(0.25, 4), rep(1, 3))

# import abudances
a2 <- read.table("pollinator_abundances.csv", sep=",", header=T, row.names = 1 )
a <- as.matrix(a2)
a[which(a>0)] <- 1

p <- read.table("plant_abundances_bin.csv", sep=",", header=T, row.names = 1 )


# arranging interactions into a list adjacency matrices (one for each site)
read.table("interactions.csv", sep=",", header=T, row.names = 1 ) %>%
  split(., .$Site) %>%
  lapply(., function(x){
    net <- matrix(x$Links, nrow=length(unique(x$Pol_sp)), ncol=length(unique(x$Plant_sp)))
    colnames(net)<- unique(x$Plant_sp)
    rownames(net) <- unique(x$Pol_sp)
    return(net)
    }) -> ntw



# calculating FD outputs
# ######################
# 
pol.coords <- species.coords(t3, a, weight)
pol.measures <- FD_measures(pol.coords$coords, pol.coords$centr, a)

pl.coords <- species.coords(pl_t, p, pl.weight)
pl.measures <- FD_measures(pl.coords$centr, pl.coords$centr, p)


# with weighted centroid: only originality should change
w.pol.coords <- species.coords(t3,a2,weight)
w.pol.measures <- FD_measures(w.pol.coords$coords, w.pol.coords$centr, a2)
colnames(w.pol.measures) <- c("w.orig", "w.redund", "w.uniq")




# Calculating bipartite metrics
###############################
ntw %>% lapply(., getspe, measure="hs") %>% do.call(c,.) -> pol.hs
ntw <- lapply(ntw, t)
ntw %>% lapply(., getspe, measure="hs") %>% do.call(c,.) -> pl.hs
ntw %>% 
  lapply(., specieslevel, index=c(
    "normalised degree", "species strength","species specificity", "d", "PDI"), 
    level="higher") %>% 
  do.call(rbind, .) -> splvl.pollis

ntw %>% 
  lapply(., specieslevel, index=c(
    "normalised degree", "species strength","species specificity", "d", "PDI"), 
    level="lower") %>% 
  do.call(rbind, .) -> splvl.plants


# removing the lassor outlier from spslevel but keeping line for later match with POLLIS 
splvl.pollis[4,] <- NA
pol.hs[4] <- NA
pol.measures <- rbind(pol.measures[1:3,], c(NA, NA, NA), pol.measures[4:89,])
w.pol.measures <- rbind(w.pol.measures[1:3,], c(NA, NA, NA), w.pol.measures[4:89,])
# combine and save data
pol.rev <- cbind(pol.measures, w.pol.measures, splvl.pollis, pol.hs)
plant.rev <- cbind(pl.measures, splvl.plants, pl.hs)
# save(pol.rev, plant.rev, file="pol_plant_FD_bip_n_insects.RData")

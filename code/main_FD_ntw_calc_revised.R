#################################################################################
# code by Camille Coux
# last modified 20 Jan 2016
#
# FD and bipartite metric calculation with ntw weighed by the number of insect 
# individuals from a same species carrying pollen of 1 plant. 
#################################################################################

library(ESM)
library(bipartite)
library(magrittr)

source("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/code/nearest_neighbors_func.R")

### 1. Import data:
### ###############

## pollinator traits
t3 <- read.table("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/data/pollinator_traits.csv", sep=",", header=T, row.names = 1 )

# assign pollinator trait weights
weight <- c(0.5, 0.5,1,0.5,0.5,1, rep(1/6, 6), rep(1,3))

## plant traits
pl_t <-  read.table("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/data/plant_traits.csv", sep=",", header=T, row.names = 1 )

# assign plant trait weights
pl.weight <- c(rep(1, 8), rep(0.25, 4), rep(1, 3))

## import abudances
a2 <- read.table("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/data/pollinator_abundances.csv", sep=",", header=T, row.names = 1 )
a <- as.matrix(a2)
a[which(a>0)] <- 1

p <- read.table("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/data/plant_abundances_bin.csv", sep=",", header=T, row.names = 1 )

# outlier check (methods from the R Book, Crawley M.J., 2007, p363)
leverage<-function(x){1/length(x)+(x-mean(x))^2/sum((x-mean(x))^2)}
x1 <- a2[a2>0]
# pdf("C:/Users/camille/Desktop/PhD/Australian data/mes trucs/FD/fdisp/proper intnat/scaled_strength/revisions/figs/final_figs/S2.pdf")
plot(leverage(x1),type="h", ylab="Pollinator abundances")
points(leverage(x1))


# A general rule is that a point is highly influential if its leverage is 
# greater than 2P/N, where p is the number of parameters in the model and 
# N is the number of data points.
P <- 5; N <- length(x1)
abline(((2*P)/N),0,lty=2)
# dev.off()
# point 47 is clearly an outlier; it corresponds to the L_sordidum abundance in
# the first site. We replace its abundance by 0 in the abundance matrix, such that 
# the FD space is calculated without it. The analysis is repeated without changing 
# the 
# abundance of this outlier and presented in Appendix 3.

# a2[1,"L_sordidum"] <- 0


## arranging interactions into a list adjacency matrices (one for each site)
interactions <- read.table("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/data/interactions.csv", sep=",", header=T, row.names = 1 )
read.table("C:/Users/camille/Desktop/Network_functional_Roles_manuscript/data/interactions.csv", sep=",", header=T, row.names = 1 ) %>%
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
pol.coords <- species.coords(t3, a, weight, na.rm=T) # species scores in the trait space
pol.measures <- FD_measures(pol.coords$coords, pol.coords$centr, a2) # keep weighted
# abundance matrix to keep track of abundances in the analysis

pl.coords <- species.coords(pl_t, p, pl.weight)
pl.measures <- FD_measures(pl.coords$coords, pl.coords$centr, p)


# with weighted centroid: only originality should change
w.pol.coords <- species.coords(t3,a2,weight)
w.pol.measures <- FD_measures(w.pol.coords$coords, w.pol.coords$centr, a2)
colnames(w.pol.measures)[which(colnames(w.pol.measures)=="orig")] <- "w.orig" 
colnames(w.pol.measures)[which(colnames(w.pol.measures)=="uniq")] <- "w.uniq" 




# Calculating bipartite metrics
###############################

# for pollinator communities
ntw %>% lapply(., getspe, measure="hs") %>% do.call(c,.) -> pol.hs
ntw <- lapply(ntw, t)
ntw %>% lapply(., getspe, measure="hs") %>% do.call(c,.) -> pl.hs
ntw %>% 
  lapply(., specieslevel, index=c(
    "normalised degree", "species strength","species specificity", "d", "PDI"), 
    level="higher") %>% 
  do.call(rbind, .) -> splvl.pollis

# for plant communities
ntw %>% 
  lapply(., specieslevel, index=c(
    "normalised degree", "species strength","species specificity", "d", "PDI"), 
    level="lower") %>% 
  do.call(rbind, .) -> splvl.plants

# removing the lassor outlier from spslevel but keeping line in 
if (a2[1,"L_sordidum"] == 0){   
  splvl.pollis[4,] <- NA
  pol.hs[4] <- NA
  pol.measures <- rbind(pol.measures[1:3,], c(NA, NA, NA), pol.measures[4:89,])
  w.pol.measures <- rbind(w.pol.measures[1:3,], c(NA, NA, NA), w.pol.measures[4:89,])
}

# combine data
pol.rev <- cbind(pol.measures, w.pol.measures, splvl.pollis, pol.hs)
plant.rev <- cbind(pl.measures, splvl.plants, pl.hs)

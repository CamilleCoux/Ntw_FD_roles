# code by Camille Coux
# last modified 21 Dec 2015

# Analysis code.

## New runs with scaled strength, and without the abundance outlier. And no intnat.
require(nlme)
require(MuMIn)


# loading data to make unique df with weighted and unweighted fdisp together
load("/Users/camille/Desktop/PhD/Australian data/mes trucs/FD/fdisp/LME_inputs_bin.RData")
load("/Users/camille/Desktop/PhD/Australian data/mes trucs/FD/fdisp/LME_inputs.RData")
load("/Users/camille/Desktop/PhD/Australian data/mes trucs/FD/fdisp/proper intnat/scaled_strength/revisions/pol_plant_FD_bip_n_insects.RData")

POLLIS <- data.frame(farm=pollis$farm, sp=pollis$sp, pol.abun = pollis$pol.abun, 
                     pol.orig = pol.rev$orig, pol.redund = pol.rev$redund, pol.uniq = pol.rev$uniq,
                     w.pol.orig = pol.rev$w.orig, w.pol.redund = pol.rev$w.redund, w.pol.uniq = pol.rev$w.uniq,
                     normalised.degree=pol.rev$normalised.degree, species.strength=pol.rev$species.strength,
                     d=pol.rev$d, PDI = pol.rev$PDI, pol.hs = pol.rev$pol.hs)


FLOWERS <- data.frame(farm=flowers$farm, sp=flowers$sp, pl.abun = flowers$pl.abun, 
                      pl.orig = plant.rev$orig, w.pl.orig = plant.rev$w.pl.orig, 
                      pl.redund = plant.rev$redund, pl.uniq = plant.rev$uniq,
                      normalised.degree=plant.rev$normalised.degree, species.strength=plant.rev$species.strength,
                      d=plant.rev$d, PDI = plant.rev$PDI, hs = plant.rev$pl.hs)

POLLIS <- POLLIS[-4,] # outlier
sizes <- lapply(ntw, length)
POLLIS %>% split(., .$farm) -> pollis
pollis <- lapply(1:21, function(i){
  x <- pollis[[i]]
  x$size <- rep(sizes[[i]], dim(x)[1])
  return(x)
})
POLLIS <- do.call(rbind, pollis)
FLOWERS %>% split(., .$farm) -> flowers
flowers <- lapply(1:21, function(i){
  x <- flowers[[i]]
  x$size <- rep(sizes[[i]], dim(x)[1])
  return(x)
})
FLOWERS <- do.call(rbind, flowers)
# save(POLLIS, FLOWERS, file="POLLIS_FLOWERS_n_insect.RData")
rm(com_metrics, flowers, pollis, flowers_bin, pollis_bin, fdispMOD)

mumin.proc <- function(model){
  ms <- dredge(model, rank="AICc", REML=F)
  msList <- get.models(ms, subset = delta < 2, method="REML")
  if (length(msList)>1){
    output <- list(ms=ms, avg = summary(model.avg(msList)))
  }else{
    output <- ms
  }
  return(output)
}


# ORIGINALITY - ND
summary(pol.bin.orig.ND<-lme(normalised.degree ~  pol.orig* size + pol.orig*pol.abun , random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(pol.bin.orig.ND)


# weighted
summary(pol.W.orig.ND<-lme(normalised.degree ~  w.pol.orig* size + pol.abun*w.pol.orig, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(pol.W.orig.ND)


# ORIGINALITY - hs
summary(pol.bin.orig.hs<-lme( pol.hs ~  pol.orig* size + pol.orig*pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(pol.bin.orig.hs)
summary(pol.bin.orig.hs<-lme( pol.hs ~  pol.orig*pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))

# weighted
summary(pol.W.orig.hs<-lme( pol.hs ~  w.pol.orig* size + pol.abun*w.pol.orig, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(pol.W.orig.hs)
summary(pol.W.orig.hs<-lme( pol.hs ~  pol.abun*w.pol.orig, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))

# REDUNDANCY (nearest neighbour) - ND
summary(pol.redund.ND<-lme(normalised.degree ~  pol.redund* size + pol.redund*pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(pol.redund.ND)
summary(pol.redund.ND<-lme(normalised.degree ~  pol.redund + size + pol.redund, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))

# REDUNDANCY (nearest neighbour) - hs
summary(pol.redund.hs<-lme(pol.hs ~  pol.redund* size + pol.redund*pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(pol.redund.hs)
summary(pol.redund.hs<-lme(pol.hs ~  pol.redund*pol.abun + size, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))



# just wondering about this since the beginning...
summary(red.worig.nd<-lme(normalised.degree ~  (pol.redund +w.pol.orig+size+pol.abun)^2, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(red.worig.nd)
summary(red.worig.nd<-lme(normalised.degree ~  pol.redund *w.pol.orig+size+pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))

summary(red.worig.nd<-lme(normalised.degree ~  (pol.redund +pol.orig+size+pol.abun)^2, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
mumin.proc(red.worig.nd)



##############
### PLANTS ###
##############

# ORIGINALITY - ND
summary(pl.bin.orig.ND<-lme(normalised.degree ~  pl.orig* size, random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))
mumin.proc(pl.bin.orig.ND)
summary(pl.bin.orig.ND<-lme(normalised.degree ~  size + pl.abun , random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))

# ORIGINALITY - HS
summary(pl.bin.orig.hs<-lme(hs ~  pl.orig* size, random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))
mumin.proc(pl.bin.orig.hs)


# REDUNDANCY (nearest neighbour) - ND
summary(pl.redund.ND<-lme(normalised.degree ~  pl.redund* size, random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))
mumin.proc(pl.redund.ND)
# nothing

# REDUNDANCY (nearest neighbour) - hs
summary(pl.redund.hs<-lme(hs ~  pl.redund* size, random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))
mumin.proc(pl.redund.hs)
# nothing

# just wondering about this since the beginning...
summary(pl.red.worig.nd<-lme(normalised.degree ~  (pl.redund +pl.orig+size)^2, random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))
mumin.proc(pl.red.worig.nd)



##############################################################################
# JACKNIFE REMOVAL OF TRAITS WITH FLIPPED MODELS.
##############################################################################

# plant traits, polli traits. original models. Remove traits one by one. recalc
# originality and uniqueness. recalc model each time and compare fit with 
# original. assign the AIC differences to traits and save them in a table. 
# should i only do this where there are significant relationships in the 
# original model? at least start with those...

# other problem lies in the model approximation... I'm gonna take the closest
# possible model conatining all variables from the averaged results.
# 
# loading the FD metics functions
source("/Users/camille/Desktop/PhD/Australian data/mes trucs/FD/fdisp/proper intnat/scaled_strength/revisions/nearest_neighbors_func.R")


# loading the trait data
load("/Users/camille/Desktop/PhD/Australian data/mes trucs/FD/fdisp/fdispMOD_inputs.RData")

### reorganising the traits table for pollis

t3 <- t2[,-3] # getting rid of body depth as super high correlation with weight
t<-t2[,14:24]
season <- NA
for (i in 1:16){season <- c(season, colnames(t)[which(t[i,] == max(t[i,]))])}
t3$season <- as.factor(season[-c(1, 4,19)]) # taking out NA and ties.
t <- t2[,26:33]
daily <- NA
for (i in 1:16){daily <- c(daily, colnames(t)[which(t[i,] == max(t[i,]))])}
t3$daily <- as.factor(daily[-1])
t3 <- t3[,-c(13:23, 25:32)] # removing superfluous cols

rm(t, daily, season, traits.sel)
weight <- c(0.5, 0.5,1,0.5,0.5,1, rep(1/6, 6), rep(1,3))


### editing plant trait table: can't get peak seasonal abund since binary.
### just removing nat/int col
pl_t <- pl_t[,-16]
pl.weight <- pl.weight[-16]

# now organising groups according to weights.
b.size <- c(1,2); vis.time <- 3; nec.pol<-c(4,5); soc.sol <- 6; larv.diet <- c(7:12)
carrying.structre <- 13; season <- 14; daily <- 15
jack.traits <- as.data.frame(matrix(NA,4 ,8))
colnames(jack.traits) <- c("b.size", "vis.time", "nec.pol", "soc.sol", "larv.diet", "carrying.structre", "season", "daily")
count <- 1
for (trait in c(b.size, vis.time, nec.pol, soc.sol, larv.diet, carrying.structre, season, daily)) {
  t3[,-trait] %>%
    species.coords(.,a2,weight[-trait]) -> w.pol.coords.jk
  
  w.pol.measures.jk <- FD_measures(w.pol.coords.jk$coords, w.pol.coords.jk$centr, a2)
  colnames(w.pol.measures.jk) <- c("w.pol.orig.jk", "w.pol.redund.jk", "w.uniq.jk")
  POLLIS$w.pol.orig.jk <- w.pol.measures.jk$w.pol.orig.jk[-4]
  POLLIS$w.pol.redund.jk <- w.pol.measures.jk$w.pol.redund.jk[-4]
  
  pol.W.orig.ND.jk <- AIC(lme(normalised.degree ~  w.pol.orig.jk*size + pol.abun*w.pol.orig.jk, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
  pol.W.orig.hs.jk <- AIC(lme( pol.hs ~  pol.abun*w.pol.orig.jk, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
  pol.redund.ND.jk <- AIC(lme(normalised.degree ~  w.pol.redund.jk + size + w.pol.redund.jk, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
  pol.redund.hs.jk <- AIC(lme(pol.hs ~  w.pol.redund.jk*pol.abun + size, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
  
  # negative delta values indicate a worse fit for jk ==> trait contributes + to model fit when in it.
  delta.pol.W.orig.ND <- AIC(pol.W.orig.ND) - pol.W.orig.ND.jk
  delta.pol.W.orig.hs <- AIC(pol.W.orig.hs) - pol.W.orig.hs.jk
  delta.pol.redund.ND <- AIC(pol.redund.ND) - pol.redund.ND.jk
  delta.pol.redund.hs <- AIC(pol.redund.hs) - pol.redund.hs.jk
  
  jack.traits[,count] <- c(delta.pol.W.orig.ND, delta.pol.W.orig.hs, delta.pol.redund.ND, delta.pol.redund.hs)
  count <- count + 1
}

jack.traits <- jack.traits[,-c(9:12)]
write.csv(jack.traits, "jack_traits.csv")
jack.traits <- read.csv("jack_traits.csv", row.names = 1)
pdf("/Users/camille/Desktop//PhD/Australian data/mes trucs/FD/fdisp/proper intnat/scaled_strength/revisions/figs/jack_traits.pdf",
       width = 8)
par(mar=c(8, 6, 4, 2))
boxplot(jack.traits[,order(colMeans(jack.traits))], ylab="delta AIC", las=2)
abline(0,0, lty=3)
dev.off()
# OLD CODE TO KEEP IN CASE REVIEWS ASK FOR MORE

###PLANTS
# # weighted
# summary(pl.W.orig.ND<-lme(normalised.degree ~  w.pl.orig* size + pl.abun*w.pl.orig, random=~1|farm, na.action=na.fail, data=FLOWERS, method="REML"))
# mumin.proc(pl.W.orig.ND)




##### PDI #############
# # ORIGINALITY - PDI
# summary(pol.bin.orig.PDI<-lme( PDI ~  pol.orig* size + pol.orig*pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
# mumin.proc(pol.bin.orig.PDI)


# summary(pol.W.orig.PDI<-lme( PDI ~  w.pol.orig* size + pol.abun*w.pol.orig, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
# mumin.proc(pol.W.orig.PDI)


# # REDUNDANCY - PDI
# summary(pol.redund.PDI<-lme(PDI ~  pol.redund* size + pol.redund*pol.abun, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
# mumin.proc(pol.redund.PDI)
# summary(pol.redund.PDI<-lme(PDI ~  pol.redund, random=~1|farm, na.action=na.fail, data=POLLIS, method="REML"))
#######################












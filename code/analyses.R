#################################################################################
# code by Camille Coux
# last modified 20 Jan 2016
#
# This is code for the analysis ran to obtain results in our manuscript. We provide 
# it as a courtesy and in an effort of transparency and reproducibility. 
# 
# 
# NB: here, weighted abundances correspond to individuals from a same species 
# carrying pollen of 1 plant. 
#################################################################################

require(nlme)
require(MuMIn)

source("main_FD_ntw_calc_revised.R")


sizes <- lapply(ntw, length)
pol.rev %>% split(., .$farm) -> pol.rev
pol.rev <- lapply(1:21, function(i){
  x <- pol.rev[[i]]
  x$size <- rep(sizes[[i]], dim(x)[1])
  return(x)
})
pol.rev <- do.call(rbind, pol.rev)
plant.rev %>% split(., .$farm) -> plant.rev
plant.rev <- lapply(1:21, function(i){
  x <- plant.rev[[i]]
  x$size <- rep(sizes[[i]], dim(x)[1])
  return(x)
})
plant.rev <- do.call(rbind, plant.rev)

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
summary(pol.bin.orig.ND<-lme(normalised.degree ~  orig* size + orig*abundance, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))
mumin.proc(pol.bin.orig.ND)


# weighted
summary(pol.W.orig.ND<-lme(normalised.degree ~  w.orig* size + abundance*w.orig, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))
mumin.proc(pol.W.orig.ND)


# ORIGINALITY - hs
summary(pol.bin.orig.hs<-lme( pol.hs ~  orig* size + orig*abundance, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))
mumin.proc(pol.bin.orig.hs)

# weighted
summary(pol.W.orig.hs<-lme( pol.hs ~  w.orig* size + abundance*w.orig, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))
mumin.proc(pol.W.orig.hs)
summary(pol.W.orig.hs<-lme( pol.hs ~ abundance*w.orig, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))


# UNIQUENESS (nearest neighbour) - ND
summary(uniq.ND<-lme(normalised.degree ~  uniq* size + uniq*abundance, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))
mumin.proc(uniq.ND)

# UNIQUENESS (nearest neighbour) - hs
summary(uniq.hs<-lme(pol.hs ~  uniq* size + uniq*abundance, random=~1|farm, na.action=na.fail, data=pol.rev, method="REML"))
mumin.proc(uniq.hs)




##############
### PLANTS ###
##############

# ORIGINALITY - ND
summary(pl.bin.orig.ND<-lme(normalised.degree ~  pl.orig* size, random=~1|farm, na.action=na.fail, data=plant.rev, method="REML"))
mumin.proc(pl.bin.orig.ND)

# ORIGINALITY - HS
summary(pl.bin.orig.hs<-lme(hs ~  pl.orig* size, random=~1|farm, na.action=na.fail, data=plant.rev, method="REML"))
mumin.proc(pl.bin.orig.hs)


# UNIQUENESS (nearest neighbour) - ND
summary(pl.uniq.ND<-lme(normalised.degree ~  pl.uniq* size, random=~1|farm, na.action=na.fail, data=plant.rev, method="REML"))
mumin.proc(pl.uniq.ND)
# nothing

# UNIQUENESS(nearest neighbour) - hs
summary(pl.uniq.hs<-lme(hs ~  pl.uniq* size, random=~1|farm, na.action=na.fail, data=plant.rev, method="REML"))
mumin.proc(pl.uniq.hs)
# nothing



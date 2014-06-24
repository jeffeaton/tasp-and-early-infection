##########################
#### Load the library ####
##########################

dyn.load("code/rlib.so")

## Dimensions of the state space
NG <- 2; RG <- 3; DS <- 6; ART_ST <- 15;

## Wrapper functions to call functions defined in rlib.cpp
R0 <- function(theta, after.cbar.cx = FALSE, rel.priminf = NULL){
  
  out <- .C("callR0", c(as.numeric(theta), rel.priminf), after.cbar.cx, numeric(1), numeric(NG*RG*(DS-1)))
  r0 <- out[[3]]
  attr(r0, "vec") <- out[[4]]
  return(r0)
}

epiGrowthRate <- function(theta, after.cbar.cx = FALSE, rel.priminf = NULL){

  out <- .C("callEpidemicGrowthRate", c(as.numeric(theta), rel.priminf), after.cbar.cx, numeric(1), numeric(NG*RG*(DS-1)))
    
  growth.rate <- out[[3]]
  attr(growth.rate, "vec") <- out[[4]]
  return(growth.rate)
}

createNGM <- function(theta, after.cbar.cx = FALSE, rel.priminf = NULL){

  Fmat <- Vmat <- matrix(0, NG*RG*(DS-1), NG*RG*(DS-1))
  out <- .C("callCreateNGM", c(as.numeric(theta), rel.priminf), after.cbar.cx, Fmat, Vmat)
  Fmat <- t(out[[3]])   # Note: tranpose because GSL stores matrices in column order
  Vmat <- t(out[[4]])   # Note: tranpose because GSL stores matrices in column order

  E <- matrix(0, NG*(DS-1)*RG, NG*RG)
  E[1:RG, 1:RG] <- diag(3)
  E[(DS-1)*RG + 1:RG, RG + 1:RG] <- diag(3)

  ngm <- t(E) %*% Fmat %*% solve(Vmat) %*% E
  
  return(list("Fmat"=t(out[[3]]), "Vmat"=t(out[[4]]), "ngm"=ngm))
}

## Simulate calibration prevalence
sim.prev <- function(theta, rel.priminf = NULL){                     
  return(.Call("callSimPrev", c(theta, rel.priminf)))
}

## Simulate model and output full state space
sim.mod <- function(theta, outDates, rateARTinit = 0.0, stageARTelig = "none", dateART = 0.0, fracScreened = 0.0, rel.priminf = NULL){
  
  stageARTeligIdx <- switch(stageARTelig, "none" = 6, "100" = 5, "200" = 4, "350" = 3, "afterprim" = 2, "all" = 1)
  return(.Call("callSimMod", c(theta, rel.priminf), outDates, rateARTinit, stageARTeligIdx, dateART, fracScreened))
}

#########################################################
####  Functions to calculate epidemiologic outcomes  ####
#########################################################

## Note: takes output from sim.mod as argument

adult.prev <- function(sim)return(rowSums(sim$X[,,1:RG,-1,,drop=FALSE],,2)/rowSums(sim$X[,,1:RG,,,drop=FALSE],,2))
all.adult.prev <- function(sim)return(rowSums(sim$X[,,1:RG,-1,,drop=FALSE])/rowSums(sim$X[,,1:RG,,,drop=FALSE]))
rg.prev <- function(sim)return(rowSums(sim$X[,,,-1,,drop=FALSE],,3)/rowSums(sim$X,,3))

adult.inc <- function(sim)return(rowSums(sim$inc,,2)/rowSums(sim$X[,,1:RG,1,,drop=FALSE],,2))
all.adult.inc <- function(sim)return(rowSums(sim$inc)/rowSums(sim$X[,,1:RG,1,,drop=FALSE]))
rg.inc <- function(sim)return(sim$inc/rowSums(sim$X[,,1:3,1,,drop=FALSE],,3))

all.adult.annualinc <- function(sim, outDates){
  ann.idx <- floor(c(0, outDates[-length(outDates)]))
  return(tapply(rowSums(sim$inc), ann.idx, sum) / tapply(rowSums(sim$X[,,1:RG,1,,drop=FALSE]), ann.idx, mean))
}

stage.trans <- function(sim)return(apply(sim$trans[,,,2:DS,,drop=FALSE], c(1,2,4), sum)/
                                   rep(rowSums(sim$trans, ,2), times= DS-1))
all.stage.trans <- function(sim)return(apply(apply(sim$trans[,,,2:DS,,drop=FALSE], c(1,4), sum), 2, "/", rowSums(sim$trans)))
prim.trans <- function(sim)return(rowSums(sim$trans[,,,2,,drop=FALSE])/rowSums(sim$trans))

pop.prev <- function(sim)return(rowSums(sim$X[,,,-1,,drop=FALSE],,2)/rowSums(sim$X,,2))
pop.size <- function(sim)return(rowSums(sim$X))
pop.on.art <- function(sim)return(rowSums(sim$X[,,,,c(3:7, 10:14)]))
pop.art.coverage <- function(sim, elig) return(pop.on.art(sim)/pop.art.need(sim, elig))
prop.on.art <- function(sim)return(pop.on.art(sim)/pop.size(sim))
prop.art.need <- function(sim, elig)return(pop.art.need(sim)/pop.size(sim))
pop.art.need <- function(sim, elig)return(switch(elig, "350" = rowSums(sim$X[,,,4:DS,]), "200" = rowSums(sim$X[,,,5:DS,])))


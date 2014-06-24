######################
####  Load model  ####
######################

setwd("../../tasp-and-early-infection")

source("code/analysis-functions.R")
source("figures/load-posterior-distributions.R")
source("figures/sa-data.R")

######################################
####  calibration HIV prevalence  ####
######################################

outDates <- 1990:2010+0.5
fitprev.fixed1 <- lapply(fixed1.theta, sim.prev, rel.priminf=fixed1.beta)
fitprev.fixed9 <- lapply(fixed9.theta, sim.prev, rel.priminf=fixed9.beta)
fitprev.fixed26 <- lapply(fixed26.theta, sim.prev, rel.priminf=fixed26.beta)


#############################
####  no ART projection  ####
#############################

fineOutDates <- seq(1990.5, 2020.7, 0.1)

mod.noart.fixed1 <- lapply(fixed1.theta, sim.mod, fineOutDates, rel.priminf=fixed1.beta)
noart.stagetrans.fixed1 <- sapply(mod.noart.fixed1, all.stage.trans)
dim(noart.stagetrans.fixed1) <- c(length(fineOutDates), DS-1, length(fixed1.theta))
noart.stagetrans.fixed1 <- noart.stagetrans.fixed1[-1,,]
noart.inc.fixed1 <- array(10*sapply(mod.noart.fixed1, adult.inc), c(length(fineOutDates), NG, length(mod.noart.fixed1)))[seq(2, 302, 10),,]
rm(mod.noart.fixed1)

mod.noart.fixed9 <- lapply(fixed9.theta, sim.mod, fineOutDates, rel.priminf=fixed9.beta)
noart.stagetrans.fixed9 <- sapply(mod.noart.fixed9, all.stage.trans)
dim(noart.stagetrans.fixed9) <- c(length(fineOutDates), DS-1, length(fixed9.theta))
noart.stagetrans.fixed9 <- noart.stagetrans.fixed9[-1,,]
noart.inc.fixed9 <- array(10*sapply(mod.noart.fixed9, adult.inc), c(length(fineOutDates), NG, length(mod.noart.fixed9)))[seq(2, 302, 10),,]
rm(mod.noart.fixed9)

mod.noart.fixed26 <- lapply(fixed26.theta, sim.mod, fineOutDates, rel.priminf=fixed26.beta)
noart.stagetrans.fixed26 <- sapply(mod.noart.fixed26, all.stage.trans)
dim(noart.stagetrans.fixed26) <- c(length(fineOutDates), DS-1, length(fixed26.theta))
noart.stagetrans.fixed26 <- noart.stagetrans.fixed26[-1,,]
noart.inc.fixed26 <- array(10*sapply(mod.noart.fixed26, adult.inc), c(length(fineOutDates), NG, length(mod.noart.fixed26)))[seq(2, 302, 10),,]
rm(mod.noart.fixed26)


#####################
####  Figure S5  ####
#####################
library(RColorBrewer)
library(adegenet)

transp.level <- 0.2

fnPlotPrev <- function(fitprev){
  maleprev <- sapply(fitprev, function(sim) sim$prev[,1])
  femaleprev <- sapply(fitprev, function(sim) sim$prev[,2])
  gamma.unif.pr <- c(0,1)
  S2 <- 1/sum(1/sa.anc$logit.var)
  dbar <- apply((sa.anc$logit.p - logit(femaleprev[which(outDates %in% (sa.anc$year + 0.5)),]))/sa.anc$logit.var, 2, sum)*S2
  gammahat <- rnorm(length(fitprev), dbar, sqrt(S2))
  while(any(gammahat < gamma.unif.pr[1] | gammahat > gamma.unif.pr[2])) 
    gammahat <- rnorm(dim(samp.prev)[3], dbar, sqrt(S2))
  ancprev <- invlogit(logit(femaleprev) + rep(gammahat, each=length(outDates)))
  ####
  plot(outDates, rep(NA, length(outDates)), ylim=c(0, .32), main = "", xlab="", ylab="", xaxt="n", yaxt="n")
  ##
  anc.plot.dates <- outDates[outDates <= 2010.5]
  anc.plot.indices <- 1:length(anc.plot.dates)
  ##
  matlines(outDates, maleprev, col=transp("lightblue", transp.level), lty=1, lwd=0.3)
  matlines(outDates, femaleprev, col=transp("lightpink", transp.level), lty=1, lwd=0.3)
  matlines(anc.plot.dates, ancprev[which(outDates %in% anc.plot.dates),], col=transp("lightgreen", transp.level), lty=1, lwd=0.3)
  ##
  lines(outDates, rowMeans(maleprev), col="royalblue", lwd=2)
  matlines(outDates, t(apply(maleprev, 1, quantile, probs=c(.025, .975))), col="royalblue", lwd=1, lty=2)
  lines(outDates, rowMeans(femaleprev), col="deeppink", lwd=2)
  matlines(outDates, t(apply(femaleprev, 1, quantile, probs=c(.025, .975))), col="deeppink", lwd=1, lty=2)
  matlines(anc.plot.dates, t(apply(ancprev[anc.plot.indices,], 1, quantile, probs=c(.025, .975))), col="forestgreen", lwd=1, lty=3)
  ##
  with(sa.hsrc.male, points(year+.5, prevalence, col="darkblue", pch=5))
  with(sa.hsrc.male, segments(x0=year+.5, invlogit(logit.p - sqrt(logit.var)*qnorm(.975)),
                              year+.5, invlogit(logit.p + sqrt(logit.var)*qnorm(.975)),
                              col="darkblue"))
  with(sa.hsrc.female, points(year+.5, prevalence, col="deeppink4", pch=5))
  with(sa.hsrc.female, segments(x0=year+.5, invlogit(logit.p - sqrt(logit.var)*qnorm(.975)),
                                year+.5, invlogit(logit.p + sqrt(logit.var)*qnorm(.975)),
                                col="deeppink4"))
  with(sa.anc, points(year+.5, prevalence, col="darkgreen"))
  with(sa.anc, segments(x0=year+.5, invlogit(logit.p - sqrt(logit.var)*qnorm(.975)),
                        year+.5, invlogit(logit.p + sqrt(logit.var)*qnorm(.975)),
                      col="darkgreen"))
  axis(1, seq(1990.5, 2010.5, 5), seq(1990, 2010, 5), las=1)
  axis(2, 0:6*0.05, labels=FALSE, las=2)
  return(NULL)
}

fnPlotInc <- function(noart.inc){
  matplot(1990:2020+0.5, 100*noart.inc[,2,], type="l", col=transp("lightpink", transp.level), lty=1, lwd=0.3, ylim=c(0, 3.7), yaxt="n",
          las=1, ylab="HIV incidence (per 100 PY)")
  matlines(1990:2020+0.5, 100*noart.inc[,1,], col=transp("lightblue", transp.level), lty=1, lwd=0.3)
  lines(1990:2020+0.5, rowMeans(100*noart.inc[,1,]), col="royalblue", lwd=2)
  matlines(1990:2020+0.5, t(apply(100*noart.inc[,1,], 1, quantile, probs=c(.025, .975))), col="royalblue", lwd=1, lty=2)
  lines(1990:2020+0.5, rowMeans(100*noart.inc[,2,]), col="deeppink", lwd=2)
  matlines(1990:2020+0.5, t(apply(100*noart.inc[,2,], 1, quantile, probs=c(.025, .975))), col="deeppink", lwd=1, lty=2)
  axis(2, 0:7*0.5, labels=FALSE)
  return(NULL)
}

fnPlotStageTrans <- function(noart.stagetrans){
  barplot(t(rowMeans(noart.stagetrans,,2)), space=0, border=NA, xaxs="i",
          col=brewer.pal(5, "Blues")[c(5, 1:4)], xlab="", ylab="", yaxt="n", cex.axis=0.9,
          xaxs="i", mgp=c(2.2, 0.5, 0))
  axis(1, at = 0:6*50, labels=FALSE)
  text(0:3*100, -0.1, seq(1990, 2020, by=10), srt=45, xpd=NA)
  axis(2, 0:5/5, labels=FALSE)
  return(NULL)
}

pdf("figures/figureS4.pdf", h=6, w=6.5, pointsize=7)
par(mfrow=c(3,3), cex=1, mar=c(2, 1, 1, 0.5), oma=c(1, 3, 1.5, 0.5), tcl=-0.25, mgp=c(2, 0.5, 0))
####           ####
####  Panel A  ####
####           ####
fnPlotPrev(fitprev.fixed1)
mtext("No increase", cex=1.3, line = 0.5, font=2)
axis(2, 0:6*0.05, paste(0:6*5, "%", sep=""), las=2, tick=FALSE, lwd=NA)
mtext("HIV prevalence", 2, 3, las=3)
mtext("A", 3, 0.2, adj=-0.25, font=2, cex=1.5)
fnPlotPrev(fitprev.fixed9)
mtext("9.2x increase", cex=1.3, line = 0.5, font=2)
fnPlotPrev(fitprev.fixed26)
mtext("26x increase", cex=1.3, line = 0.5, font=2)
####           ####
####  Panel B  ####
####           ####
fnPlotInc(noart.inc.fixed1)
axis(2, 0:7*0.5, las=1)
mtext("HIV incidence (per 100 PY)", 2, 3, las=3)
mtext("B", 3, 0.25, adj=-0.25, font=2, cex=1.5)
fnPlotInc(noart.inc.fixed9)
fnPlotInc(noart.inc.fixed26)
####           ####
####  Panel C  ####
####           ####
fnPlotStageTrans(noart.stagetrans.fixed1)
axis(2, 0:5/5, paste(0:5*20, "%", sep=""), las=1, lwd=NA)
mtext("Percentage of transmissions", 2, 3, las=3)
mtext("C", 3, 0.25, adj=-0.25, font=2, cex=1.5)
fnPlotStageTrans(noart.stagetrans.fixed9)
fnPlotStageTrans(noart.stagetrans.fixed26)
dev.off()

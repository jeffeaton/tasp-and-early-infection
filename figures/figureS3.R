######################
####  Load model  ####
######################

setwd("../../tasp-and-early-infection")

source("code/analysis-functions.R")
source("figures/load-posterior-distributions.R")


#############################
####  no ART projection  ####
#############################

outDates <- seq(2009.9, 2040.0, 0.1)
noart.out.post <- lapply(post.theta, sim.mod, outDates)
priminf.trans.post <- sapply(noart.out.post, prim.trans)[-1,]
noart.inc.post <- sapply(noart.out.post, all.adult.inc)[-(1:2),]
rm(noart.out.post)


############################
####  ART intervention  ####
############################

sigma <- c(0, 4.166666666666667, 0.219298245614035, 0.217391304347826, 0.239644878712710, 0.973559205743437)
ART.fut <- 2010

sim.corr <- function(elig, coverage, screen.frac = 0.0){

  tx.idx <- switch(elig, "all" = 2, "afterprim" = 3, "350" = 4, "200" = 5)

  screen.cov <- coverage*screen.frac
  unscr.coverage <- 1-(1-coverage)/(1-screen.cov)
  lambda <- uniroot(function(x)return(prod(sigma[tx.idx:6]/(x+sigma[tx.idx:6])) - (1-unscr.coverage)), c(0, 10))$root

  system.time(art.out.post <- lapply(post.theta, sim.mod, outDates, lambda, elig, ART.fut, screen.frac))
  art.inc.post <- sapply(art.out.post, all.adult.inc)[-(1:2),]
  art.increduc.post <- ((noart.inc.post - art.inc.post)/noart.inc.post)
  
  return(data.frame(increduc = rowMeans(art.increduc.post), corr2010 = apply(art.increduc.post, 1, cor, priminf.trans.post[1,])))
}

cov.seq <- c(0.8, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.95)
elig.seq <- c("350", "all", "200")
scrfrac.seq <- c(0.0, 1.0)

increduc.cor.covseq <- lapply(cov.seq, sim.corr, elig="350")
increduc.cor.eligseq <- lapply(elig.seq, sim.corr, cov=0.8)
increduc.cor.scrfracseq <- lapply(scrfrac.seq, sim.corr, elig="350", cov=0.8)


#####################
####  Figure S3  ####
#####################

library(RColorBrewer)

pdf("figures/figureS3.pdf", h=4, w=6.5, pointsize=7)
par(mfrow=c(3,2), oma=c(2, 1, 1, 9), mar=c(2, 3, 1, 1), cex=1, tcl=-0.25, mgp=c(2, 0.5, 0))
##
matplot(outDates[-(1:2)] - 2010,
        1-sapply(increduc.cor.eligseq, "[[", "increduc")[,3:1],
        type="l", col=c(brewer.pal(3, "Set1")[1:2], 1),
        lty=1, lwd=2, ylim=c(0.4, 1.0), las=1, yaxt="n", ylab="")
axis(2, 4:10/10, NA, las=1)
axis(2, 2:5/5, paste(3:0*20, "%", sep=""), las=1)
mtext("% incidence reduction", 2, 2.5)
mtext("A", 3, 0.5, at=-6.5, cex=1.5, font=2, xpd=NA)
##
plot(outDates[-(1:2)] - 2010, rep(NA, 300),
     type="n", ylim=c(-1,1), las=1,
     ylab="Correlation", xlab="Years since treatment intervention",
     main="")
abline(h=c(-1, 0, 1), col="lightgrey", lty=c(1,2,1))
matlines(outDates[-(1:2)] - 2010,
         sapply(increduc.cor.eligseq, "[[", "corr2010")[,3:1],
         col=c(brewer.pal(3, "Set1")[1:2], 1),
         lwd=2, lty=1)
legend("right", legend=c("CD4 < 350", "all HIV+", "CD4 < 200"), lwd=2, lty=1, col=c(1, brewer.pal(3, "Set1")[2:1]), inset=-0.45, xpd=NA,
       title="ART eligibility")
####
matplot(outDates[-(1:2)] - 2010, 1-sapply(increduc.cor.covseq, "[[", "increduc")[, c(2:9,1)],
        type="l", col=c(brewer.pal(8, "RdBu"), 1),
        lty=1, lwd=2, ylim=c(0.4, 1.0), las=1, yaxt="n", ylab="")
axis(2, 4:10/10, NA, las=1)
axis(2, 2:5/5, paste(3:0*20, "%", sep=""), las=1)
mtext("% incidence reduction", 2, 2.5)
mtext("B", 3, 0.5, at=-6.5, cex=1.5, font=2, xpd=NA)
##
plot(outDates[-(1:2)] - 2010, rep(NA, 300),
     type="n", ylim=c(-1,1), las=1,
     ylab="Correlation", xlab="Years since treatment intervention",
     main="")
abline(h=c(-1, 0, 1), col="lightgrey", lty=c(1,2,1))
matlines(outDates[-(1:2)] - 2010,
         sapply(increduc.cor.covseq, "[[", "corr2010")[, c(2:9,1)],
         col=c(brewer.pal(8, "RdBu"), 1),
         lwd=2, lty=1)
legend("right", legend=paste(100*cov.seq, "%", sep=""), lwd=2, lty=1, col=c(1, brewer.pal(8, "RdBu")), inset=-0.5, xpd=NA, ncol=2,
       title="% accessing treatment", cex=0.9)
####
matplot(outDates[-(1:2)] - 2010, 1-sapply(increduc.cor.scrfracseq, "[[", "increduc")[,2:1],
        type="l", col=c(brewer.pal(3, "Set1")[1], 1),
        lty=1, lwd=2, ylim=c(0.4, 1.0), las=1, yaxt="n", ylab="")
axis(2, 4:10/10, NA, las=1)
axis(2, 2:5/5, paste(3:0*20, "%", sep=""), las=1)
mtext("% incidence reduction", 2, 2.5)
mtext("Years after ART intervention", 1, 2)
mtext("C", 3, 0.5, at=-6.5, cex=1.5, font=2, xpd=NA)
##
plot(outDates[-(1:2)] - 2010, rep(NA, 300),
     type="n", ylim=c(-1,1), las=1,
     ylab="Correlation", xlab="Years since treatment intervention",
     main="")
abline(h=c(-1, 0, 1), col="lightgrey", lty=c(1,2,1))
matlines(outDates[-(1:2)] - 2010, sapply(increduc.cor.scrfracseq, "[[", "corr2010")[,2:1],
         col=c(brewer.pal(3, "Set1")[1], 1),
         lwd=2, lty=1)
legend("right", legend=c("constant rate\n", "mean 1 year\nafter eligibility"), lwd=2, lty=1, col=c(1, brewer.pal(3, "Set1")[1]), inset=-0.48, xpd=NA,
       title="Timing of initiation")
mtext("Years after ART intervention", 1, 2)
dev.off()

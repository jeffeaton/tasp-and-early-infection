######################
####  Load model  ####
######################

setwd("../../tasp-and-early-infection")

source("code/analysis-functions.R")
source("figures/load-posterior-distributions.R")
source("figures/sa-data.R")


############################################
####  model calibration HIV prevalence  ####
############################################

outDates <- 1990:2010+0.5
fitprev.post <- lapply(post.theta, sim.prev)

male.prev <- sapply(fitprev.post, function(sim) sim$prev[,1])
female.prev <- sapply(fitprev.post, function(sim) sim$prev[,2])

gamma.unif.pr <- c(0,1)
S2 <- 1/sum(1/sa.anc$logit.var)
dbar <- apply((sa.anc$logit.p - logit(female.prev[which(outDates %in% (sa.anc$year + 0.5)),]))/sa.anc$logit.var, 2, sum)*S2
gamma.hat <- rnorm(length(post.theta), dbar, sqrt(S2))
while(any(gamma.hat < gamma.unif.pr[1] | gamma.hat > gamma.unif.pr[2])) ## this is inefficient (braindead even) if gamma is  close to prior limits...
  gamma.hat <- rnorm(length(post.theta), dbar, sqrt(S2))
anc.prev <- invlogit(logit(female.prev) + rep(gamma.hat, each=length(outDates)))


#############################
####  no ART projection  ####
#############################

fineOutDates <- seq(1990.5, 2020.7, 0.1)
mod.noart.post <- lapply(post.theta, sim.mod, fineOutDates)
no.art.inc <- array(10*sapply(mod.noart.post, adult.inc), c(length(fineOutDates), NG, length(post.theta)))[seq(2, 302, 10),,]
noart.combinc.post <- sapply(mod.noart.post, all.adult.inc)[-1,]

system.time(noart.stagetrans.post <- sapply(mod.noart.post, all.stage.trans))
dim(noart.stagetrans.post) <- c(length(fineOutDates), DS-1, length(post.theta))
noart.stagetrans.post <- noart.stagetrans.post[-1,,]
priminf.trans2010.post <- noart.stagetrans.post[202, 1,]


####################
####  Figure 1  ####
####################
library(RColorBrewer)
library(adegenet)

transp.level <- 0.2

pdf("figures/figure1.pdf", h=(5+3.8)*0.3937, w=8.7*0.3937, pointsize=7)
par(oma=c(0,0,0,0))
layout(rbind(1, rep(2:3, each=3)), h=c(5, 3.8))
par(mar=c(3, 4.5, 0.5, 1.5), mgp=c(2, 0.5, 0), tcl=-0.25, cex=1)
####           ####
####  Panel A  ####
####           ####
outDates <- 1990:2010+0.5
plot(outDates, rep(NA, length(outDates)), ylim=c(0, .32), main = "", xlab="", ylab="", xaxt="n", yaxt="n")
##
anc.plot.dates <- outDates[outDates <= 2010.5]
anc.plot.indices <- 1:length(anc.plot.dates)
##
matlines(outDates, male.prev, col=transp("lightblue", transp.level), lty=1, lwd=0.3)
matlines(outDates, female.prev, col=transp("lightpink", transp.level), lty=1, lwd=0.3)
matlines(anc.plot.dates, anc.prev[which(outDates %in% anc.plot.dates),], col=transp("lightgreen", transp.level), lty=1, lwd=0.3)
##
lines(outDates, rowMeans(male.prev), col="royalblue", lwd=1.5)
matlines(outDates, t(apply(male.prev, 1, quantile, probs=c(.025, .975))), col="royalblue", lwd=1, lty=2)
lines(outDates, rowMeans(female.prev), col="deeppink", lwd=1.5)
matlines(outDates, t(apply(female.prev, 1, quantile, probs=c(.025, .975))), col="deeppink", lwd=1, lty=2)
matlines(anc.plot.dates, t(apply(anc.prev[anc.plot.indices,], 1, quantile, probs=c(.025, .975))), col="forestgreen", lwd=1, lty=3)
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
mtext("HIV prevalence", 2, 2.5, las=3)
axis(2, 0:6*0.05, paste(0:6*5, "%", sep=""), las=2)
axis(1, seq(1990.5, 2010.5, 5), seq(1990, 2010, 5), las=1)
legend("bottomright", legend=c("posterior mean", "95% credible interval"),
       lty=c(1,2), inset=c(0, 0.02), lwd=c(1.5,1), bty="n", title="Posterior prevalence", cex=0.9)
legend("topleft", legend=c("ANC Prevalence", "Females, HSRC Survey", "Males, HSRC Survey"),
         pch=c(1,5,5), col=c("darkgreen", "deeppink4", "darkblue"), lty=1, inset=c(0, 0.02),
         title="Data mean and 95% CI", bty="n", cex=0.8)
mtext("A", 3, -1.2, adj=-0.18, font=2, cex=1.5)
####           ####
####  Panel B  ####
####           ####
par(mar=c(1.5,3,1.0,0.5), cex=1)
matplot(1990:2020+0.5, 100*no.art.inc[,2,], type="l", col=transp("lightpink", transp.level), lty=1, lwd=0.3,
        las=1, ylab="HIV incidence (per 100 PY)", xlab="")
matlines(1990:2020+0.5, 100*no.art.inc[,1,], col=transp("lightblue", transp.level), lty=1, lwd=0.3)
lines(1990:2020+0.5, rowMeans(100*no.art.inc[,1,]), col="royalblue", lwd=1.5)
matlines(1990:2020+0.5, t(apply(100*no.art.inc[,1,], 1, quantile, probs=c(.025, .975))), col="royalblue", lwd=1, lty=2)
lines(1990:2020+0.5, rowMeans(100*no.art.inc[,2,]), col="deeppink", lwd=1.5)
matlines(1990:2020+0.5, t(apply(100*no.art.inc[,2,], 1, quantile, probs=c(.025, .975))), col="deeppink", lwd=1, lty=2)
legend("topright", c("Females", "Males", "mean", "95% CI"), pch=c(15, 15, NA, NA), pt.cex=c(1.5, 1.5, NA, NA), lwd=c(NA, NA, 1.5, 1), lty=c(NA, NA, 1, 2),
       col=c("deeppink", "royalblue", 1, 1), cex=0.7, bty="n")
mtext("B", 3, 0.5, adj=-0.26, font=2, cex=1.5)
####           ####
####  Panel C  ####
####           ####
cd4.stages <- c("Early Inf.", "CD4 >350", "CD4 200-350", "CD4 100-200", "CD4 <100")
par(mar=c(1.5,3,1.0,0.8), cex=1)
barplot(t(rowMeans(noart.stagetrans.post,,2)), space=0, border=NA, xaxs="i",
        col=brewer.pal(5, "Blues")[c(5, 1:4)], xlab="Year", ylab="Percentage of transmissions", yaxt="n", cex.axis=0.9,
        xaxs="i", mgp=c(2.2, 0.5, 0))
axis(1, at = 0:6*50, labels=FALSE)
axis(1, at = 0:3*100, labels=seq(1990, 2020, by=10), tick=FALSE, cex.axis=0.9)
axis(2, 0:5/5, paste(0:5*20, "%", sep=""), cex.axis=0.9, las=1)
legend("bottomleft", legend=cd4.stages[5:1], pch=15, col=brewer.pal(5, "Blues")[c(4:1, 5)], bg="white", inset=c(.02, .01), cex=.6, pt.cex=1.5)
mtext("C", 3, 0.5, adj=-0.28, font=2, cex=1.5)
dev.off()


####                   ####
###  Analysis for text  ###
####                   ####

## peak incidence
mean(fineOutDates[apply(noart.combinc.post, 2, which.max) + 1])

## reduction over the next decade
mean(1 - noart.combinc.post[apply(noart.combinc.post, 2, which.max) + 100 + (0:(ncol(post.theta) - 1))*nrow(noart.combinc.post)] / apply(noart.combinc.post, 2, max))
quantile(1 - noart.combinc.post[apply(noart.combinc.post, 2, which.max) + 100 + (0:(ncol(post.theta) - 1))*nrow(noart.combinc.post)] / apply(noart.combinc.post, 2, max), c(0.025, 0.975))

## reduction in contact rate
mean(as.numeric(post.theta[14,]))
quantile(as.numeric(post.theta[14,]), c(0.025, 0.975))


####################
####  Figure 2  ####
####################

library(adegenet)

pdf("figures/figure2.pdf", h=4*0.3937, w=8.7*0.3937, pointsize=7)
par(mfrow=c(1,2), mar=c(3, 3, 1.4, 1), mgp=c(1.7, 0.5, 0), tcl=-0.25, cex=1)
####           ####
####  Panel A  ####
####           ####
plot(density(priminf.trans2010.post), type="n", ylim=c(0, 15), 
     main="", xlab="Early transmission in 2010", ylab="Density", xaxt="n",
     las=1, xpd=FALSE)
abline(h=0, col="grey", lwd=0.5)
lines(density(priminf.trans2010.post), lwd=1.5, col=3)
lines(density(r0current.post, from=1.1, to=8.0), lwd=1.5, col=4, xlab=c(0, 5))
axis(1, 1:6/20, labels=FALSE, tick=TRUE, cex.axis=0.9)
axis(1, 1:3/10, paste(1:3*10, "%", sep=""), tick=FALSE, cex.axis=0.9)
mtext("A", 3, -0.1, adj=-0.28, font=2, cex=1.5)
####           ####
####  Panel B  ####
####           ####
plot(log(as.numeric(post.theta[17,])), priminf.trans2010.post, pch=16, cex=0.2, col=transp("gray30", 0.4), ylim=c(0.05, 0.3),
     xaxt="n", yaxt="n", xlab="Relative early infect. (log scale)", ylab="")
axis(1, log(c(7, 10, 15, 20, 30, 50)), paste(c(7, 10, 15, 20, 30, 50), "x", sep=""), cex.axis=0.9)
axis(2, 1:6/20, paste(1:6*5, "%", sep=""), tick=TRUE, cex.axis=0.9, las=1)
mtext("Early transmission", 2, 2.2)
abline(lm(priminf.trans2010.post~log(as.numeric(post.theta[17,]))), col=2)
mtext(bquote(R^2 == .(round(summary(lm(priminf.trans2010.post~log(as.numeric(post.theta[17,]))))$r.squared, 2))), 3, -1.3, adj=0.05, font=2, cex=1)
mtext("B", 3, -0.1, adj=-0.32, font=2, cex=1.5)
dev.off()


####                   ####
###  Analysis for text  ###
####                   ####

## early transmission
mean(priminf.trans2010.post)
quantile(priminf.trans2010.post, c(0.025, 0.975))

post.theta <- as.matrix(post.theta)
summary(lm(priminf.trans2010.post ~ log(post.theta[17,])))$r.squared
summary(lm(priminf.trans2010.post ~ log(post.theta[17,]) + post.theta[13,]))$r.squared

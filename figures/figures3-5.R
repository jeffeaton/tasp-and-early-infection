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
noart.annualinc.post <- sapply(noart.out.post, all.adult.annualinc, outDates)[-(1:2),]
rm(noart.out.post)


### fixed early infectiousness projections (Figure 5B and Table 1 only)

noart.out.fixed1 <- lapply(fixed1.theta, sim.mod, outDates, rel.priminf=fixed1.beta)
noart.inc.fixed1 <- sapply(noart.out.fixed1, all.adult.inc)
priminf.trans.fixed1 <- sapply(noart.out.fixed1, prim.trans)[-1,]
rm(noart.out.fixed1)

noart.out.fixed9 <- lapply(fixed9.theta, sim.mod, outDates, rel.priminf=fixed9.beta)
noart.inc.fixed9 <- sapply(noart.out.fixed9, all.adult.inc)
priminf.trans.fixed9 <- sapply(noart.out.fixed9, prim.trans)[-1,]
rm(noart.out.fixed9)

noart.out.fixed26 <- lapply(fixed26.theta, sim.mod, outDates, rel.priminf=fixed26.beta)
noart.inc.fixed26 <- sapply(noart.out.fixed26, all.adult.inc)
priminf.trans.fixed26 <- sapply(noart.out.fixed26, prim.trans)[-1,]
rm(noart.out.fixed26)


############################
####  ART intervention  ####
############################

sigma <- c(0, 4.166666666666667, 0.219298245614035, 0.217391304347826, 0.239644878712710, 0.973559205743437)
ART.fut <- 2010
coverage <- 0.8
lam350 <- uniroot(function(x)return(prod(sigma[4:6]/(x+sigma[4:6])) - (1-coverage)), c(0, 10))$root

art350.out.post <- lapply(post.theta, sim.mod, outDates, lam350, "350", ART.fut)
art350.inc.post <- sapply(art350.out.post, all.adult.inc)[-(1:2),]
art350.annualinc.post <- sapply(art350.out.post, all.adult.annualinc, outDates)[-(1:2),]
art350.increduc.post <- ((noart.inc.post - art350.inc.post)/noart.inc.post)
art350.annualincreduc.post <- ((noart.annualinc.post - art350.annualinc.post)/noart.annualinc.post)
rm(art350.out.post)


### fixed early infectiousness projections (Figure 5B and Table 1 only)

art350.out.fixed1 <- lapply(fixed1.theta, sim.mod, outDates, lam350, "350", ART.fut, rel.priminf=fixed1.beta)
art350.inc.fixed1 <- sapply(art350.out.fixed1, all.adult.inc)
art350.increduc.fixed1 <- ((noart.inc.fixed1 - art350.inc.fixed1)/noart.inc.fixed1)
rm(art350.out.fixed1)

art350.out.fixed9 <- lapply(fixed9.theta, sim.mod, outDates, lam350, "350", ART.fut, rel.priminf=fixed9.beta)
art350.inc.fixed9 <- sapply(art350.out.fixed9, all.adult.inc)
art350.increduc.fixed9 <- ((noart.inc.fixed9 - art350.inc.fixed9)/noart.inc.fixed9)
rm(art350.out.fixed9)

art350.out.fixed26 <- lapply(fixed26.theta, sim.mod, outDates, lam350, "350", ART.fut, rel.priminf=fixed26.beta)
art350.inc.fixed26 <- sapply(art350.out.fixed26, all.adult.inc)
art350.increduc.fixed26 <- ((noart.inc.fixed26 - art350.inc.fixed26)/noart.inc.fixed26)
rm(art350.out.fixed26)


##########################
####  R0 calculation  ####
##########################

r0start.post <- sapply(post.theta, R0)
r0current.post <- sapply(post.theta, R0, TRUE)

### fixed early infectiousness projections (Figure 5B and Table 1 only)

r0start.fixed1 <- sapply(fixed1.theta, R0, rel.priminf=fixed1.beta)
r0current.fixed1 <- sapply(fixed1.theta, R0, TRUE, rel.priminf=fixed1.beta)
r0start.fixed9 <- sapply(fixed9.theta, R0, rel.priminf=fixed9.beta)
r0current.fixed9 <- sapply(fixed9.theta, R0, TRUE, rel.priminf=fixed9.beta)
r0start.fixed26 <- sapply(fixed26.theta, R0, rel.priminf=fixed26.beta)
r0current.fixed26 <- sapply(fixed26.theta, R0, TRUE, rel.priminf=fixed26.beta)


####################
####  Figure 3  ####
####################

library(adegenet) # for transparent colors [transp()]

psmp <- 1:length(post.theta)  # in case we want to down sample to plot points
priminf.xlim <- c(0.1, 0.26)
transp.level <- 0.4

pdf("figures/figure3.pdf", h=9*0.3937, w=8.7*0.3937, pointsize=7)
par(oma=c(0,0,0.25,0))
layout(rbind(1:2, 3:4, 5), h=c(1,1,1.2))
par(tcl=-0.25, cex=1.0, cex.axis=0.8)
####           ####
####  Panel A  ####
####           ####
par(mar=c(1, 2+3.25, 0, 0.25), mgp=c(2, 0.4, 0), tcl=-0.25, cex=1.0)
## 
plot(priminf.trans.post[1,psmp], art350.annualincreduc.post[1,psmp], pch=16, cex=0.25, col=transp("gray40", transp.level), xlim=priminf.xlim, yaxt="n", xaxt="n", ylab="")
mtext("Year 1", 3, -1.15, adj=0.95, font=2, cex=0.9)
axis(2, seq(0.055, 0.07, 0.0025), labels=FALSE, cex.axis=0.8)
axis(2, seq(0.055, 0.07, 0.005), sprintf("%.1f%%", 100*seq(0.055, 0.07, 0.005)), tick=FALSE, las=1, cex.axis=0.8)
axis(1, seq(0.05, 0.3, 0.05), NA, las=1, cex.axis=0.8)
abline(lm(art350.annualincreduc.post[1,]~priminf.trans.post[1,]), col="red")
mtext("% reduction in incidence", 2, 2.5, xpd=NA, adj=8.5)
mtext("% transmission during early infection", 3, -21, outer=TRUE, cex=1, padj=-1)
mtext("A", 3, -1.1, at=0.005, cex=1.5, font=2)
##
par(mar=c(1, 2, 0, 3.25))
plot(priminf.trans.post[1,psmp], art350.annualincreduc.post[2,psmp], pch=16, cex=0.25, col=transp("gray40", transp.level), xlim=priminf.xlim, ylim=c(0.145, 0.19), yaxt="n", xaxt="n")
mtext("Year 2", 3, -1.15, adj=0.95, font=2, cex=0.9)
axis(2, seq(0.12, 0.19, 0.01), labels=FALSE, tick=TRUE, las=1, cex.axis=0.8)
axis(2, seq(0.13, 0.19, 0.01), paste(100*seq(0.13, 0.19, 0.01), "%", sep=""), tick=FALSE, las=1, cex.axis=0.8)
axis(1, seq(0.05, 0.3, 0.05), NA, las=1, cex.axis=0.8)
abline(lm(art350.annualincreduc.post[2,]~priminf.trans.post[1,]), col="red")
##
par(mar=c(1, 2+3.25, 0, 0.25))
plot(priminf.trans.post[1,psmp], art350.annualincreduc.post[10,psmp], pch=16, cex=0.25, col=transp("gray40", transp.level), xlim=priminf.xlim, yaxt="n", xaxt="n", ylab="")
mtext("Year 10", 3, -1.15, adj=0.95, font=2, cex=0.9)
axis(2, seq(0.15, 0.4, 0.025), labels=FALSE, tick=TRUE, las=1, cex.axis=0.8)
axis(2, seq(0.15, 0.4, 0.05), paste(100*seq(0.15, 0.4, 0.05), "%", sep=""), tick=FALSE, las=1, cex.axis=0.8)
axis(1, seq(0.05, 0.3, 0.05), paste(100*seq(0.05, 0.3, 0.05), "%", sep=""), las=1, cex.axis=0.8)
abline(lm(art350.annualincreduc.post[10,]~priminf.trans.post[1,]), col="red")
##
par(mar=c(1, 2, 0, 3.25))
plot(priminf.trans.post[1,psmp], art350.annualincreduc.post[30,psmp], pch=16, cex=0.25, col=transp("gray40", transp.level), xlim=priminf.xlim, ylim=c(0.1, 0.35), yaxt="n", xaxt="n")
mtext("Year 30", 3, -1.15, adj=0.95, font=2, cex=0.9)
axis(2, seq(0.1, 0.5, 0.05), labels=FALSE, tick=TRUE, las=1, cex.axis=0.8)
axis(2, seq(0.1, 0.5, 0.05), paste(100*seq(0.1, 0.5, 0.05), "%", sep=""), tick=FALSE, las=1, cex.axis=0.8)
axis(1, seq(0.05, 0.3, 0.05), paste(100*seq(0.05, 0.3, 0.05), "%", sep=""), las=1, cex.axis=0.8)
abline(lm(art350.annualincreduc.post[30,]~priminf.trans.post[1,]), col="red")
####           ####
####  Panel B  ####
####           ####
par(mar=c(2.5, 3.5, 2.5, 1.5), cex.axis=0.8, mgp=c(1.3, 0.4, 0))
plot(seq(0.1, 30, 0.1), apply(art350.increduc.post, 1, cor, priminf.trans.post[1,]),
     type="n", ylim=c(-1.0,1.0), las=1, lwd=2,
     ylab="Correlation", xlab="Years after treatment intervention",
     main="", xaxt="n", xpd=NA)
axis(1, seq(0, 30, 5))
abline(h=c(-1, 0, 1), col="lightgrey", lty=c(1,2,1))
lines(seq(0.1, 30, 0.1), apply(art350.increduc.post, 1, cor, priminf.trans.post[1,]),
      type="l", las=1, lwd=2)
mtext("B", 3, -0, at=-5, cex=1.5, font=2)
dev.off()


####                   ####
###  Analysis for text  ###
####                   ####


rowMeans(art350.annualincreduc.post[c(1, 2, 10, 30),])
apply(art350.annualincreduc.post[c(1, 2, 10, 30),], 1, quantile, c(0.025, 0.975))
apply(art350.annualincreduc.post[c(1, 2, 10, 30),], 1, cor, priminf.trans.post[1,])



####################
####  Figure 4  ####
####################

pdf("figures/figure4.pdf", h=8*0.3937, w=8.7*0.3937, pointsize=7)
par(oma=c(0, 0, 0.5, 0))
layout(rbind(c(1, 1, 2, 2), c(0, 3, 3, 0)))
par(mar=c(2.7, 3, 0.9, 1), mgp=c(1.7, 0.5, 0), tcl=-0.25, cex=1)
####           ####
####  Panel A  ####
####           ####
plot(density(r0start.post, from=1.1, to=8.0), lwd=1.5, col=3, xlim=c(1.1, 8.0), ylim=c(0, 0.8),
     main="", xlab=expression(R[0]), ylab="",
     las=1, xpd=FALSE)
lines(density(r0current.post, from=1.1, to=8.0), lwd=1.5, col=4, xlab=c(0, 5))
mtext("Density", 2, 2)
legend("topleft", c("current  ", "start"), col=4:3, pch=15, pt.cex=2, horiz=TRUE, cex=0.75, x.intersp=1, bty="n")
mtext("A", 3, 0, adj=-0.3, font=2, cex=1.5)
####           ####
####  Panel B  ####
####           ####
plot(r0start.post, priminf.trans.post[1,], pch=16, cex=0.2, col=transp("green4", 0.3),
     xlim=c(3, 8), ylim=c(0.08, 0.32),
     las=1, yaxt="n", ylab="", xlab="", xpd=NA)
axis(2, 1:6/20, FALSE, cex.axis=0.9, las=1)
axis(2, 1:6/20, paste(1:6*5, "%", sep=""), lwd=NA, cex.axis=0.9, las=1)
mtext(expression(R[0]~"(start)"), 1, 1.7)
mtext("Early transmission", 2, 2.2)
mtext(bquote(R^2 == .(round(summary(lm(priminf.trans.post[1,]~I(1/r0start.post)))$r.squared, 2))), 3, -1.3, adj=0.95, font=2, cex=1)
lines(seq(1, 8, 0.1), c(lm(priminf.trans.post[1,]~I(1/r0start.post))$coef %*% rbind(1, 1/seq(1, 8, 0.1))), col=2)
mtext("B", 3, 0, adj=-0.3, font=2, cex=1.5)
####           ####
####  Panel C  ####
####           ####
par(mar=c(2.7, 2, 0.9, 2))
plot(r0current.post, art350.increduc.post[300,], pch=16, cex=0.2, col=transp("blue4", 0.3),
     xlim=c(2, 6), ylim=c(0.1, 0.35), yaxt="n", xaxt="n", ylab="", xlab=expression(R[0]~"(current)"), cex.axis=0.9)
axis(2, 0:5/10, paste(0:5*10, "%", sep=""), tick=FALSE, las=1, cex.axis=0.9)
axis(2, 0:10/20, labels=FALSE, las=1, cex.axis=0.9)
axis(1, 2:6, 2:6, cex.axis=0.9)
mtext(bquote(R^2 == .(round(summary((lm(art350.increduc.post[300,]~I(1/r0current.post))))$r.squared, 2))), 3, -1.3, adj=0.95, font=2, cex=1)
lines(seq(1, 8, 0.1), c(lm(art350.increduc.post[300,]~I(1/r0current.post))$coef %*% rbind(1, 1/seq(1, 8, 0.1))), col=2)
mtext("Reduction in incidence\n after 30 years", 2, 2.0)
mtext("C", 3, -0.2, adj=-0.4, font=2, cex=1.5)
dev.off()


####                   ####
###  Analysis for text  ###
####                   ####

mean(r0start.post); quantile(r0start.post, c(0.025, 0.975))
mean(r0current.post); quantile(r0current.post, c(0.025, 0.975))

summary(lm(priminf.trans.post[1,]~I(1/r0start.post)))$r.squared

summary(lm(r0start.post ~ t(log(post.theta[17,])) + t(post.theta[13,])))$r.squared
summary(lm(r0current.post ~ t(log(post.theta[17,])) + t(post.theta[13,])))$r.squared
summary(lm(r0current.post ~ t(log(post.theta[17,])) + t(post.theta[13,]) + t(post.theta[14,])))$r.squared

summary((lm(art350.increduc.post[300,]~I(1/r0current.post))))$r.squared



####################
####  Figure 5  ####
####################

library(RColorBrewer)
library(adegenet)

psmp <- 1:length(fixed1.theta)  # in case we want to down sample plot points
transp.level <- 0.4

pdf("figures/figure5.pdf", h=8*0.3937, w=8.7*0.3937, pointsize=7)
par(mfrow=c(2,1), mgp=c(2.7, 0.5, 0), tcl=-0.25, mar=c(3.2, 5, 1.0, 2))
####           ####
####  Panel A  ####
####           ####
plot(log(as.numeric(post.theta[17,])), art350.increduc.post[300,],
     pch=16, cex=0.25, col=transp("gray20", transp.level), ylim=c(0.07, 0.4), yaxt="n", xaxt="n",
     xlab="", ylab="")
abline(lm(art350.increduc.post[300,]~log(as.numeric(post.theta[17,]))), col="red")
mtext(bquote(R^2 == .(sprintf("%.2f", summary(lm(art350.increduc.post[300,]~log(as.numeric(post.theta[17,]))))$r.squared))), 3, -1.3, adj=0.98, font=2, cex=1)
axis(1, log(c(7, 10, 15, 20, 30, 50)), paste(c(7, 10, 15, 20, 30, 50), "x", sep=""))
axis(2, 0:5/10, paste(0:5*10, "%", sep=""), las=1)
mtext("Relative infectiousness during early infection (log scale)", 1, 1.5)
mtext("Reduction in HIV incidence rate after 30 years", 2, -1.1, outer=TRUE)
mtext("A", 3, -0.2, at = 1.44, cex=1.5, font=2, xpd=NA)
####           ####
####  Panel B  ####
####           ####
par(mar=c(2.7, 5, 0.75, 2))
set.seed(69696076)
plot(jitter(rep(1:3, each=length(psmp)), 1.5),
     c(art350.increduc.fixed1[300,psmp],
       art350.increduc.fixed9[300,psmp],
       art350.increduc.fixed26[300,psmp]),
     pch=16, cex=.25, col=transp(rep(brewer.pal(4, "Set1")[c(3,1,2)], each=length(psmp)), transp.level),
     xaxt="n", yaxt="n", main="",
     xlab="", ylab="",
     ylim=c(0.07, 0.4), xlim=c(0.6, 3.4))
segments(y0=mean(art350.increduc.fixed1[300,]), x0=0.6, x1=1.4, col="darkgreen", lwd=1.5)
segments(y0=mean(art350.increduc.fixed9[300,]), x0=1.6, x1=2.4, col="darkred", lwd=1.5)
segments(y0=mean(art350.increduc.fixed26[300,]), x0=2.6, x1=3.4, col="darkblue", lwd=1.5)
axis(2, at=0:6/10, paste(0:6*10, "%", sep=""), las=1)
axis(1, 1:3, c("No increase", "9.2x increase", "26x increase"), padj=.5, tcl=-0.1, cex.axis=0.95, mgp=c(2, 0, 0))
mtext("Relative infectiousness during early infection", 1, 1.5)
mtext("B", 3, -0.2, at = 0.08, cex=1.5, font=2, xpd=NA)
dev.off()


###################
####  Table 1  ####
###################

print.interval <- function(x, percentage=TRUE, num.dig=1, probs=c(0.025, 0.975)){
  if(percentage)
    sprintf("%.1f%% (%.1f--%.1f%%)", mean(x), quantile(x, probs[1]), quantile(x, probs[2]))
  else
    sprintf("%.1f (%.1f--%.1f)", mean(x), quantile(x, probs[1]), quantile(x, probs[2]))
}

cbind("early.trans.2010"=c(print.interval(100*priminf.trans.fixed1[1,]),
        print.interval(100*priminf.trans.fixed9[1,]),
        print.interval(100*priminf.trans.fixed26[1,])),
      "r0.start"=c(print.interval(r0start.fixed1, FALSE),
        print.interval(r0start.fixed9, FALSE),
        print.interval(r0start.fixed26, FALSE)),
      "r0.current"=c(print.interval(r0current.fixed1, FALSE),
        print.interval(r0current.fixed9, FALSE),
        print.interval(r0current.fixed26, FALSE)),
      "inc.reduc.2040"=c(print.interval(100*art350.increduc.fixed1[300,]),
        print.interval(100*art350.increduc.fixed9[300,]),
        print.interval(100*art350.increduc.fixed26[300,])))

####                   ####
###  Analysis for text  ###
####                   ####

mean(art350.increduc.post[300,])

apply(art350.increduc.post, 1, cor, log(as.numeric(post.theta[17,])))

sort(abs(cor(r0start.post, t(post.theta))))
cor(r0current.post, t(post.theta))

c(mean(art350.increduc.fixed1[300,]), quantile(art350.increduc.fixed1[300,], c(0.025, 0.975)))
c(mean(art350.increduc.fixed9[300,]), quantile(art350.increduc.fixed9[300,], c(0.025, 0.975)))
c(mean(art350.increduc.fixed26[300,]), quantile(art350.increduc.fixed26[300,], c(0.025, 0.975)))

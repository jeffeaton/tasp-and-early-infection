######################
####  Load model  ####
######################

setwd("../../tasp-and-early-infection")

source("code/analysis-functions.R")
source("figures/load-posterior-distributions.R")


#############################
####  no ART projection  ####
#############################

outDates <- seq(2009.9, 2010.0, 0.1)
noart.out.post <- lapply(post.theta, sim.mod, outDates)
priminf.trans.post <- sapply(noart.out.post, prim.trans)[-1,]
rm(noart.out.post)

r0start.post <- sapply(post.theta, R0)
r0current.post <- sapply(post.theta, R0, TRUE)

priminf <- as.numeric(post.theta[17,])
psi <- as.numeric(post.theta[13,])
creduc <- as.numeric(post.theta[14,])
  
post.theta <- as.matrix(post.theta)


#####################
####  Figure S2  ####
#####################
library(adegenet)
transp.level <- 0.3

pdf("figures/figureS2.pdf", h=5.5, w=6.5, pointsize=7)
par(oma=c(2.5, 0, 0, 0))
layout(rbind(c(1:2, 0), c(3:4, 0), 5:7))
## layout.show(7)
par(cex=1, cex.axis=0.9, mar=c(0.5, 4, 1.5, 0.5), tcl=-0.25, mgp=c(2, 0.5, 0))
####           ####
####  Panel A  ####
####           ####
plot(log(priminf), priminf.trans.post, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(0.05, 0.3),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(priminf.trans.post~log(priminf)), col=2)
mtext(bquote(R^2 == .(round(summary(lm(priminf.trans.post~log(priminf)))$r.squared, 2))), 3, -1.3, adj=0.05, font=2, cex=1)
axis(1, log(c(7, 10, 15, 20, 30, 50)), paste(c(7, 10, 15, 20, 30, 50), "x", sep=""), cex.axis=0.9)
axis(2, 1:6/20, paste(1:6*5, "%", sep=""), tick=TRUE, cex.axis=0.9, las=1)
mtext("Early transmission in 2010", 2, 2.5)
mtext("A", 3, 0, adj=-0.27, font=2, cex=1.5)
##
plot(psi, lm(priminf.trans.post~log(priminf))$res, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(-0.06, 0.08),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(lm(priminf.trans.post~log(priminf))$res ~ psi), col=2)
mtext(bquote("Residual R"^2 == .(round(summary(lm(lm(priminf.trans.post~log(priminf))$res ~ psi))$r.squared, 2))), 3, -1.3, adj=0.05, font=2, cex=1)
mtext(bquote("Combined R"^2 == .(round(summary(lm(priminf.trans.post~log(priminf)+psi))$r.squared, 2))), 3, -2.3, adj=0.05, font=
  2, cex=1)
axis(2, -3:4/50, paste(-3:4*2, "%", sep=""), tick=TRUE, cex.axis=0.9, las=1)
axis(1, cex.axis=0.9)
mtext("Residual of early transmission\nadj. for early infectiousness", 2, 2)
####           ####
####  Panel B  ####
####           ####
plot(log(priminf), r0start.post, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(3, 8),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(r0start.post~log(priminf)), col=2)
mtext(bquote(R^2 == .(round(summary(lm(r0start.post~log(priminf)))$r.squared, 2))), 3, -1.3, adj=0.95, font=2, cex=1)
axis(1, log(c(7, 10, 15, 20, 30, 50)), paste(c(7, 10, 15, 20, 30, 50), "x", sep=""), cex.axis=0.9)
axis(2, cex.axis=0.9, las=1)
mtext(expression(R[0]~"(start)"), 2, 2)
mtext("B", 3, 0, adj=-0.27, font=2, cex=1.5)
##
plot(psi, lm(r0start.post~log(priminf))$res, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(-2, 2),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(lm(r0start.post~log(priminf))$res ~ psi), col=2)
mtext(bquote("Residual R"^2 == .(round(summary(lm(lm(r0start.post~log(priminf))$res ~ psi))$r.squared, 2))), 3, -1.3, adj=0.95, font=2, cex=1)
mtext(bquote("Combined R"^2 == .(round(summary(lm(r0start.post~log(priminf)+psi))$r.squared, 2))), 3, -2.3, adj=0.95, font=
  2, cex=1)
axis(2, cex.axis=0.9, las=1)
axis(1, cex.axis=0.9)
mtext(expression("Residual of"~R[0]~"(start)"),2, 2.5)
mtext("adj. for early infectiousness", 2, 1.65)
####           ####
####  Panel C  ####
####           ####
plot(log(priminf), r0current.post, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(2, 7),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(r0current.post~log(priminf)), col=2)
mtext(bquote(R^2 == .(round(summary(lm(r0current.post~log(priminf)))$r.squared, 2))), 3, -1.3, adj=0.95, font=2, cex=1)
axis(1, log(c(7, 10, 15, 20, 30, 50)), paste(c(7, 10, 15, 20, 30, 50), "x", sep=""), cex.axis=0.9)
axis(2, cex.axis=0.9, las=1)
mtext(expression(R[0]~"(current)"), 2, 2)
mtext("C", 3, 0, adj=-0.27, font=2, cex=1.5)
mtext("Relative early infectiousness ("~beta[E]/beta[A]~")", 1, 2, cex=1)
##
plot(psi, lm(r0current.post~log(priminf))$res, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(-2, 2.5),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(lm(r0current.post~log(priminf))$res ~ psi), col=2)
mtext(bquote("Residual R"^2 == .(round(summary(lm(lm(r0current.post~log(priminf))$res ~ psi))$r.squared, 2))), 3, -1.3, adj=0.95, font=2, cex=1)
mtext(bquote("Combined R"^2 == .(round(summary(lm(r0current.post~log(priminf)+psi))$r.squared, 2))), 3, -2.3, adj=0.95, font=
  2, cex=1)
axis(2, cex.axis=0.9, las=1)
axis(1, cex.axis=0.9)
mtext(expression("Residual of"~R[0]~"(current)"),2, 2.5)
mtext("adj. for early infectiousness", 2, 1.65)
mtext(expression("Rate move high to low risk ("~psi~")"), 1, 2, cex=1)
##
plot(creduc, lm(r0current.post~log(priminf))$res, pch=16, cex=0.3, col=transp("gray30", transp.level), ylim=c(-2, 2.5),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(lm(lm(r0current.post~log(priminf))$res ~ creduc), col=2)
mtext(bquote("Residual R"^2 == .(sprintf("%.2f", summary(lm(lm(r0current.post~log(priminf))$res ~ creduc))$r.squared))), 3, -1.3, adj=0.95, font=2, cex=1)
mtext(bquote("Combined R"^2 == .(round(summary(lm(r0current.post~log(priminf)+creduc))$r.squared, 2))), 3, -2.3, adj=0.95, font=
  2, cex=1)
axis(2, cex.axis=0.9, las=1)
axis(1, 0:5/10, paste(0:5*10, "%", sep=""), cex.axis=0.9)
mtext(expression("Residual of"~R[0]~"(current)"),2, 2.5)
mtext("adj. for early infectiousness", 2, 1.65)
mtext(expression("Reduction in contact rate ("~Delta[bar(c)]~")"), 1, 2, cex=1)
dev.off()

## legend text
summary(lm(r0current.post~log(priminf)+psi+creduc))$r.squared

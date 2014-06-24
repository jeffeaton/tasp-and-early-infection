######################
####  Load model  ####
######################

setwd("../../tasp-and-early-infection")

source("code/analysis-functions.R")
source("figures/load-posterior-distributions.R")
source("figures/sa-data.R")

## Simulate the prevalence for calibration
outDates <- 1990:2010+0.5
fitprev.post <- lapply(post.theta, sim.prev)
maleprev.post <- sapply(fitprev.post, function(sim) sim$prev[,1])
femaleprev.post <- sapply(fitprev.post, function(sim) sim$prev[,2])

## Sample from posterior distribution of gamma
gamma.unif.pr <- c(0,1)
S2 <- 1/sum(1/sa.anc$logit.var)
dbar <- apply((sa.anc$logit.p - logit(femaleprev.post[which(outDates %in% (sa.anc$year + 0.5)),]))/sa.anc$logit.var, 2, sum)*S2
gammahat.post <- rnorm(length(post.theta), dbar, sqrt(S2))
while(any(gammahat.post < gamma.unif.pr[1] | gammahat.post > gamma.unif.pr[2])) ## this is inefficient (braindead even) if gamma is close to prior limits...
  gammahat.post <- rnorm(length(post.theta), dbar, sqrt(S2))

## Sample from the prior distribution
prior.samp <- cbind(sample.prior(1e7), runif(1e7, gamma.unif.pr[1], gamma.unif.pr[2]))


#####################
####  Figure S1  ####
#####################

param.idx <- c(1, 2, 3, 4, 5, 13, 6, 14, 15, 16, 7, 8, 9, 10, 11, 12, 17, 18)
param.names <- c(expression(t[0]), expression(1- pi[L]^M), expression(1- pi[L]^F),
                 expression(pi[H]^M/(1-pi[L]^M)), expression(pi[H]^F/(1-pi[L]^F)), expression(bar(c)),
                 expression(omega[M]^F), expression(omega[H]^F - omega[M]^F), expression(epsilon),
                 expression(kappa[H]), expression(kappa[M]), expression(kappa[L]),
                 expression(psi), expression(Delta[bar(c)]), expression(t[bar(c)]), expression(t[bar(c)] + delta[bar(c)]),
                 expression(beta[E]/beta[A]), expression(gamma))
param.label.adj <- c(0.05, 0.95, 0.05, 0.05, 0.95, 0.05, 0.05, 0.95, 0.05, 0.95, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.95, 0.05)

post.theta.mat <- as.matrix(rbind(post.theta, gammahat.post))
pr.limits <- rbind(theta.unif.pr, gamma.unif.pr)

pdf("figures/figureS1.pdf", h=6.5, w=6.5, pointsize=8)
layout(t(matrix(c(1:19, 19), 4, 5)))
par(las=1, mgp=c(2, 0.5, 0), tcl=-0.3, mar=c(2,2.2,0.5,0.7), oma=c(0,0,0,0), cex=1)
for(i in param.idx){
  plot(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2]),
       xlim=pr.limits[i,], ylab="", xlab="", main="",
       ylim=1.2*c(0, max(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y)))
  lines(density(prior.samp[,i], from = pr.limits[i,1], to = pr.limits[i,2]), col=2, lty=2, lwd=0.5)
  mtext(param.names[i],3,-1.4, adj=param.label.adj[i], font=2, cex=1.0)
  segments(x0 = median(post.theta.mat[i,]), y0=0,
           y1=max(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y)*1.05,
           col="green3", lwd=1.5)
  segments(x0 = quantile(post.theta.mat[i,], c(0.025, 0.975)), y0=0,
           y1=max(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y)*1.05,
           col="blue", lwd=1, lty=2)
}
par(mar=c(0, 0, 0, 0), cex=1)
plot(0, 0, type="n", bty="n", axes=FALSE)
legend("center", c("Posterior density", "Prior density", "Posterior mean", "95% CI"), col=c(1, 2, "green3", "blue"), lty=c(1, 2, 1, 1), lwd=c(1, 0.5, 1.5, 1), bty="n", cex=1.2)
dev.off()

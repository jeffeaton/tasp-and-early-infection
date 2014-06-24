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
fitprev.post <- lapply(post.theta, sim.prev)
fitprev.fixed1 <- lapply(fixed1.theta, sim.prev, rel.priminf=fixed1.beta)
fitprev.fixed9 <- lapply(fixed9.theta, sim.prev, rel.priminf=fixed9.beta)
fitprev.fixed26 <- lapply(fixed26.theta, sim.prev, rel.priminf=fixed26.beta)

sample.gamma <- function(fitprev){
  maleprev <- sapply(fitprev, function(sim) sim$prev[,1])
  femaleprev <- sapply(fitprev, function(sim) sim$prev[,2])
  S2 <- 1/sum(1/sa.anc$logit.var)
  dbar <- apply((sa.anc$logit.p - logit(femaleprev[which(outDates %in% (sa.anc$year + 0.5)),]))/sa.anc$logit.var, 2, sum)*S2
  gammahat <- rnorm(length(fitprev), dbar, sqrt(S2))
  while(any(gammahat < gamma.unif.pr[1] | gammahat > gamma.unif.pr[2])) ## this is inefficient (braindead even) if gamma is close to the prior limits...
    gammahat <- rnorm(length(fitprev), dbar, sqrt(S2))
  return(gammahat)
}

gammahat.post <- sample.gamma(fitprev.post)
gammahat.fixed1 <- sample.gamma(fitprev.fixed1)
gammahat.fixed9 <- sample.gamma(fitprev.fixed9)
gammahat.fixed26 <- sample.gamma(fitprev.fixed26)

### sample from the prior distribution
prior.samp <- cbind(sample.prior(1e7), runif(1e7, gamma.unif.pr[1], gamma.unif.pr[2]))


#####################
####  Figure S5  ####
#####################

param.idx <- c(1, 2, 3, 4, 5, 13, 6, 14, 15, 16, 7, 8, 9, 10, 11, 12, 17, 18)
param.names <- c(expression(t[0]), expression(1- pi[L]^M), expression(1- pi[L]^F),
                 expression(pi[H]^M/(1-pi[L]^M)), expression(pi[H]^F/(1-pi[L]^F)), expression(bar(c)),
                 expression(omega[M]^F), expression(omega[H]^F - omega[M]^F), expression(epsilon),
                 expression(kappa[H]), expression(kappa[M]), expression(kappa[L]),
                 expression(psi), expression(Delta[bar(c)]), expression(t[bar(c)]), expression(t[bar(c)] + delta[bar(c)]),
                 expression(beta[E]/beta[A]), expression(gamma))
param.label.adj <- c(0.05, 0.95, 0.05, 0.05, 0.95, 0.05, 0.05, 0.95, 0.05, 0.95, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.95, 0.05)

pr.limits <- rbind(theta.unif.pr, gamma.unif.pr)
post.theta.mat <- as.matrix(rbind(post.theta, gammahat.post))
fixed1.theta.mat <- as.matrix(rbind(fixed1.theta, 0, gammahat.fixed1))
fixed9.theta.mat <- as.matrix(rbind(fixed9.theta, 0, gammahat.fixed9))
fixed26.theta.mat <- as.matrix(rbind(fixed26.theta, 0, gammahat.fixed26))

pdf("figures/figureS5.pdf", h=6.5, w=6.5, pointsize=8)
layout(t(matrix(c(1:19, 19), 4, 5)))
par(las=1, mgp=c(2, 0.5, 0), tcl=-0.3, mar=c(2,2.2,0.5,0.7), oma=c(0,0,0,0), cex=1)
for(i in param.idx){
  plot(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2]),
       xlim=pr.limits[i,], ylab="", xlab="", main="", col="grey60", lwd=0.5,
       ylim=c(0, 1.2*max(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y,
         density(fixed1.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y,
         density(fixed9.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y,
         density(fixed26.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y)))
  lines(density(prior.samp[,i], from = pr.limits[i,1], to = pr.limits[i,2]), col=2, lty=2, lwd=0.5)
  if(i != 17){
    lines(density(fixed1.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2]), col="darkgreen")
    lines(density(fixed9.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2]), col="darkred")
    lines(density(fixed26.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2]), col="darkblue")
  } else {
    max.yval <- max(density(post.theta.mat[i,], from = pr.limits[i,1], to = pr.limits[i,2])$y)*1.05
    segments(x0 = fixed1.beta, y0=0, y1=max.yval, col="darkgreen", lwd=1.5)
    segments(x0 = fixed9.beta, y0=0, y1=max.yval, col="darkred", lwd=1.5)
    segments(x0 = fixed26.beta, y0=0, y1=max.yval, col="darkblue", lwd=1.5)
  } 
  mtext(param.names[i],3,-1.4, adj=param.label.adj[i], font=2, cex=0.9)
}
par(mar=c(0, 0, 0, 0), cex=1)
plot(0, 0, type="n", bty="n", axes=FALSE)
legend("center", c("1x increased infectiousness", "9.2x increased infectiousness", "26 times increased infectiousness", "Full posterior (Fig. S1)", "Prior density"),
       col=c("darkgreen", "darkred", "darkblue", "grey60", "red"), lty=c(1, 1, 1, 1, 2), lwd=c(1, 1, 1, 0.5, 0.5), bty="n", cex=1.2)
dev.off()

####################
####  Table S3  ####
####################

round(cor(t(rbind(post.theta, gammahat.post))[,param.idx]), 2)
t(t(sapply(cor(t(rbind(post.theta, gammahat.post))[,param.idx]), function(vec) sprintf("\\textcolor[gray]{%f}{%.1f}", 0.95*(1-abs(vec)), vec))))

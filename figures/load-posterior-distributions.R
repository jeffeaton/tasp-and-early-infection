set.seed(718614023); post.theta <- data.frame(t(read.table("output/postResample.txt", colClasses="numeric")[sample(1e5, 1e4),]))
set.seed(833677571); fixed1.theta <- data.frame(t(read.table("output/fixed1Resample.txt", colClasses="numeric")[sample(1e5, 1e4),]))
set.seed(1610965797); fixed9.theta <- data.frame(t(read.table("output/fixed9Resample.txt", colClasses="numeric")[sample(1e5, 1e4),]))
set.seed(396525859); fixed26.theta <- data.frame(t(read.table("output/fixed26Resample.txt", colClasses="numeric")[sample(1e5, 1e4),]))

fixed1.beta <- 1.0
fixed9.beta <- 9.2
fixed26.beta <- 26.0
  
theta.unif.pr <- cbind(c(1983.0, 0.05, 0.05, 0.2, 0.2, 0.5,  1.0,  0.0, 0.2, 0.0, 0.0, 0.0, 0.00, 0.0, 1990.0, 2002.0, 1),
                       c(1988.0, 0.70, 0.70, 0.8, 0.8, 4.0, 70.0, 50.0, 0.8, 1.0, 1.0, 1.0, 0.15, 0.7, 2002.0, 2010.0, 60))
priminf.rel.pr.mu <- 3.2
priminf.rel.pr.sigma <- 0.34
gamma.unif.pr <- c(0,1)

sample.prior <- function(n){
  mat <- t(array(runif(n*nrow(theta.unif.pr), theta.unif.pr[,1], theta.unif.pr[,2]), c(nrow(theta.unif.pr), n)))
  mat[,12] <- mat[,12]^(1/3)
  mat[,11] <- mat[,11]^(1/2)*mat[,12]
  mat[,10] <- mat[,10]*mat[,11]
  mat[,17] <- rlnorm(n, priminf.rel.pr.mu, priminf.rel.pr.sigma)
  return(mat)
}

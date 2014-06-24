theta <- c(1985.1294593287, 0.2063071211, 0.3767595068, 0.4118538427, 0.2459994280, 1.7867638957, 29.3724160181, 20.7896279814, 0.4413274774, 0.2403828905, 0.4245210215, 0.7609544656, 0.0087297793, 0.1883205148, 1997.6939461648, 2010.6496928662)

r0 <- 0

source("analysis-functions.R")

R0(theta, FALSE)
R0(theta, TRUE)

.Call("callSimMod", theta, outDates))

test <- .Call("callSimMod", theta, 0.0, 4, 2010.5, 1990:2040 + 0.5)

out <- sim.mod(theta, 1990:2040 + 0.5)
fit <- sim.prev(theta)

adult.prev(out)
all.stage.trans(out)
prop.on.art(out)

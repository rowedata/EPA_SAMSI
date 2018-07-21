
source("R/functions.R")

##Reading AQS Grid Data
AQS.grid <- read.csv("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data/AQS_DATA.csv")
NW.NR <- read.csv("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data/combinedOfNW_NR.csv")

NW.NR.AQS <- matrix(c(AQS.grid$yi, AQS.grid$yi), ncol = 2)
NW.NR.AQS.mu <- matrix(c(NW.NR$m1, NW.NR$m2), ncol = 2)
NW.NR.AQS.sd <- matrix(c(NW.NR$s1, NW.NR$s2), ncol = 2)
NW.NR.AQS.d <- matrix(c(NW.NR$d.to.R_1.edge, NW.NR$d.to.R_2.edge), ncol = 2)

##Compute MLE
phi.mle <- (coef(mle(likelihoodFunDensity, start = list(phi = 1), fixed = list(
  Y = NW.NR.AQS, dist = NW.NR.AQS.d, mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd))))[1]

##Read all DS estimates for region
NW.NR.DS <- read.csv("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data/combinedOverlap.csv")

##Calculate distance from boundaries 
DS.d.bound <- matrix(c(118 + NW.NR.DS$Longitude, -1*NW.NR.DS$Longitude - 109), ncol = 2)

##Calculate distance from centre
DS.d.cent <- abs(113.5 + NW.NR.DS$Longitude)

##Read DS estimates and standard errors
DS.y <- matrix(c(NW.NR.DS$Prediction, NW.NR.DS$Prediction.1), ncol = 2)
DS.se <- matrix(c(NW.NR.DS$SEpred, NW.NR.DS$SEpred.1), ncol = 2)

##Get estimates based on MLE
spliced.NW.NR.D <- smoothEstimate(DS.d.bound, DS.y, phi.mle)

##Plot combined density function
x <- seq(-50, 100, length = 1000)
NW.pdf <- dnorm(x, mean = NW.NR.DS$Prediction[2], sd = NW.NR.DS$SEpred[2])
NR.pdf <- dnorm(x, mean = NW.NR.DS$Prediction.1[2], sd = NW.NR.DS$SEpred.1[2])
w.1 <- exp(-1*phi.mle*DS.d.bound[2, 1])
w.2 <- exp(-1*phi.mle*DS.d.bound[2, 2])
spliced.pdf <- (w.1/(w.1 + w.2))*NW.pdf + w.2/(w.1 + w.2)*NR.pdf
par(mfrow = c(3, 1))
plot(x, NW.pdf, type = "l", ylab = "NW")
plot(x, NR.pdf, type = "l", ylab = "NR")
plot(x, spliced.pdf, type = "l", ylab = "Spliced")

##Get smoothed estimates from mixture model
phi.mle.2 <- (coef(mle(combineRVLL, start = list(phi = .4), fixed = list(
  Y = NW.NR.AQS2, dist = NW.NR.AQS.d, mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd))))[1]

spliced.NW.NR.RV <- smoothEstimate(DS.d.bound, DS.y, phi.mle2)

##Get MLE for phi which linearly varies with distance from centre
estimates.3 <- (coef(mle(likelihoodFunRV2, start = list(a1 = 1, a2 = 1), fixed = list(
  Y = NW.NR.AQS, dist = NW.NR.AQS.d, mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd, dist.cen = DS.d.cent))))


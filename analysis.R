
source("R/functions.R")

AQS.grid <- read.csv("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data/AQS_DATA.csv")
NW.NR <- read.csv("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data/combinedOfNW_NR.csv")

NW.NR.AQS <- matrix(c(AQS.grid$yi, AQS.grid$yi), ncol = 2)
NW.NR.AQS.mu <- matrix(c(NW.NR$m1, NW.NR$m2), ncol = 2)
NW.NR.AQS.sd <- matrix(c(NW.NR$s1, NW.NR$s2), ncol = 2)
NW.NR.AQS.d <- matrix(c(NW.NR$d.to.R_1.edge, NW.NR$d.to.R_2.edge), ncol = 2)

phi.mle <- (coef(mle(likelihoodFunction, start = list(phi = 1), fixed = list(
  Y = NW.NR.AQS, dist = NW.NR.AQS.d, mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd))))[1]

NW.NR.DS <- read.csv("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data/combinedOverlap.csv")

DS.d <- matrix(c(118 + NW.NR.DS$Longitude, -1*NW.NR.DS$Longitude - 109), ncol = 2)
DS.y <- matrix(c(NW.NR.DS$Prediction, NW.NR.DS$Prediction.1), ncol = 2)

spliced.NW.NR <- smoothEstimate(DS.d, DS.y, phi.mle)

x <- seq(-50, 100, length = 1000)
NW.pdf <- dnorm(x, mean = NW.NR.DS$Prediction[2], sd = NW.NR.DS$SEpred[2])
NR.pdf <- dnorm(x, mean = NW.NR.DS$Prediction.1[2], sd = NW.NR.DS$SEpred.1[2])
w.1 <- exp(-1*phi.mle*DS.d[2, 1])
w.2 <- exp(-1*phi.mle*DS.d[2, 2])
spliced.pdf <- (w.1/(w.1 + w.2))*NW.pdf + w.2/(w.1 + w.2)*NR.pdf
par(mfrow = c(3, 1))
plot(x, NW.pdf, type = "l", ylab = "NW")
plot(x, NR.pdf, type = "l", ylab = "NR")
plot(x, spliced.pdf, type = "l", ylab = "Spliced")


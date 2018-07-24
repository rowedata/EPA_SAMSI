
source("R/functions.R")

##Reading AQS Grid Data
AQS.grid <- read.csv("data/AQS_DATA.csv")
NW.NR <- read.csv("data/combinedOfNW_NR.csv")

NW.NR.AQS <- matrix(c(AQS.grid$yi, AQS.grid$yi), ncol = 2)
NW.NR.AQS2 <- c(AQS.grid$yi)
NW.NR.AQS.mu <- matrix(c(NW.NR$m1, NW.NR$m2), ncol = 2)
NW.NR.AQS.sd <- matrix(c(NW.NR$s1, NW.NR$s2), ncol = 2)
NW.NR.AQS.d <- matrix(c(NW.NR$d.to.R_1.edge, NW.NR$d.to.R_2.edge), ncol = 2)

##Compute MLE
phi.mle <- (coef(mle(likelihoodFunDensity, start = list(phi = 1), fixed = list(
  Y = NW.NR.AQS, dist = NW.NR.AQS.d, mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd))))[1]

##Read all DS estimates for region
NW.NR.DS <- read.csv("data/combinedOverlap.csv")

##Calculate distance from boundaries 
DS.d.bound <- matrix(c(118 + NW.NR.DS$Longitude, -1*NW.NR.DS$Longitude - 109)
                     , ncol = 2)

##Calculate distance from centre
DS.d.cent.lon <- abs(113.5 + NW.NR.DS$Longitude)
DS.d.cent.lat <- abs(45.5 - NW.NR.DS$Latitude)

#Centre of NW, NR
NW.centre <- c((51 + 40)/2, -(127 + 109)/2) #45.5, 118
NR.centre <- c((51 + 38)/2, -(118 + 93)/2) #44.5, 105.5

#Latitudinal and Longitudinal distance from centre of NW, NR for AQS sites
AQS.NW.distance <- matrix(c(abs(NW.centre[1] - AQS.grid$Latitude), 
                     abs(NW.centre[2] - AQS.grid$Longitude)), ncol = 2)

AQS.NR.distance <- matrix(c(abs(NR.centre[1] - AQS.grid$Latitude), 
                        abs(NR.centre[2] - AQS.grid$Longitude)), ncol = 2)

NW.NR.lat.centre.dist <- matrix(c(AQS.NW.distance[, 1], AQS.NR.distance[, 1])
                                , ncol = 2)
NW.NR.lon.centre.dist <- matrix(c(AQS.NW.distance[, 2], AQS.NR.distance[, 2])
                                , ncol = 2)

#Latitudinal and Longitudinal distance from centre of NW, NR for DS estimates
DS.NW.distance <- matrix(c(abs(NW.centre[1] - NW.NR.DS$Latitude), 
                           abs(NW.centre[2] - NW.NR.DS$Longitude)), ncol = 2)

DS.NR.distance <- matrix(c(abs(NR.centre[1] - NW.NR.DS$Latitude), 
                            abs(NR.centre[2] - NW.NR.DS$Longitude)), ncol = 2)

NW.NR.DS.lat.centre <- matrix(c(DS.NW.distance[, 1], DS.NR.distance[, 1]), ncol = 2)
NW.NR.DS.lon.centre <- matrix(c(DS.NW.distance[, 2], DS.NR.distance[, 2]), ncol = 2)

##calculate distance from center for AQS
match.site <- match(AQS.grid$Loc_Label1, NW.NR.DS$Loc_Label1)

NW.NR.AQS.cent <- DS.d.cent.lon[match.site]

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
estimates.3 <- coef(mle(likelihoodFunRV2, start = list(a1 = 0.1, a2 = 0.2), 
                         fixed = list(Y = NW.NR.AQS2, dist = NW.NR.AQS.d, 
                          mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd, 
                          dist.cen = NW.NR.AQS.cent)))

spliced.NW.NR.RV.new <- smoothEstimate2(a1 = estimates.3[1], a2 = estimates.3[2],
                                        dist.cen = DS.d.cent.lon,
                                        dist = DS.d.bound, mu = DS.y)

NW.NR.DS$spliced <- spliced.NW.NR.RV.new

estimates.4 <- coef(mle(likelihoodFunRV3, start = list(beta.0 = 0, beta.1 = 0,
                                                       alpha.0 = 0, alpha.1 = 0),
                        fixed = list(d.1.centre = NW.NR.AQS.cent[, 1], 
                                     d.2.centre = NW.NR.AQS.cent[, 2], 
                                     d.1 = NW.NR.lat.centre.dist,
                                     d.2 = NW.NR.lon.centre.dist,
                                     mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd,
                                     Y = NW.NR.AQS2)))

spliced.NW.NR.RV.4 <- smoothEstimate3(beta.0 = estimates.4[1], 
                                      beta.1 = estimates.4[2],
                                      alpha.0 = estimates.4[3],
                                      alpha.1 = estimates.4[4],
                                      d.1.centre = DS.d.cent.lat,
                                      d.2.centre = DS.d.cent.lat,
                                      d.1 = NW.NR.DS.lat.centre,
                                      d.2 = NW.NR.DS.lon.centre,
                                      mu = DS.y)


###################################
##North West and West

NW.W.AQS <- read.csv("data/NW_W.csv")
NW.W.DS <- read.csv("data/DSOverlap_NW_W.csv")

#Get DS estimates and standard errors
NW.W.AQS.mu <- NW.W.AQS[, c(6, 8)]
NW.W.AQS.sd <- NW.W.AQS[, c(7, 9)]

#Get distance from centre of intersection (lat = 42)
NW.W.AQS.d.cent <- abs(NW.W.AQS$Latitude - 42)
NW.W.DS.d.cent <- abs(NW.W.DS$Latitude - 42)

#Get all DS estimates 
NW.W.DS.mu <- NW.W.DS[, c(5, 7)]

#Get distances from boundaries of intersection (-127, -112, 44, 40)
NW.W.AQS.d.bound <- matrix(c(abs(40 - NW.W.AQS$Latitude), 
                             abs(44 - NW.W.AQS$Latitude)), ncol = 2)

NW.W.DS.d.bound <- matrix(c(abs(40 - NW.W.DS$Latitude), 
                             abs(44 - NW.W.DS$Latitude)), ncol = 2)

estimates.5 <- coef(mle(likelihoodFunRV2, start = list(a1 = 0, a2 = 0),
                    fixed = list(dist.cen = NW.W.AQS.d.cent, Y = NW.W.AQS$yi, 
                                 dist = NW.W.AQS.d.bound, mu = NW.W.AQS.mu, 
                                 sd = NW.W.AQS.sd)))

spliced.NW.W <- smoothEstimate2(a1 = estimates.5[1], a2 = estimates.5[2], 
                                dist.cen = NW.W.DS.d.cent, 
                                dist = NW.W.DS.d.bound, mu = NW.W.DS.mu)

###################################
##North Rockies and West

NR.W.AQS <- read.csv("data/W_NR.csv")
NR.W.DS <- read.csv("data/DSOverlap_NR_W.csv")

#Replace NR estimates with previous run
NW.NR.DS$spliced <- spliced.NW.NR.RV.new
NR.matched <- match(NW.NR.DS$Loc_Label1, NR.W.DS$Loc_Label1)
NR.matched <- NR.matched[!is.na(NR.matched)]
NR.W.DS[NR.matched, 5] <- NW.NR.DS$spliced

for (i in 1:dim(NW.NR.DS)[1]){
  newPredictionLocation <- match(NW.NR.DS$Loc_Label1[i], NR.W.DS$Loc_Label1)
  NR.W.DS$m1[newPredictionLocation] = NW.NR.DS$spliced[i]
}

#Get DS estimates and standard errors
NR.W.AQS.mu <- NR.W.AQS[, c(8, 6)]
NR.W.AQS.sd <- NR.W.AQS[, c(9, 7)]

#Get distance from centre of intersection (lat = 41)
NR.W.AQS.d.cent <- abs(NR.W.AQS$Latitude - 41)
NR.W.DS.d.cent <- abs(NR.W.DS$Latitude - 41)

#Get all DS estimates 
NR.W.DS.mu <- NR.W.DS[, c(5, 7)]

#Get distances from boundaries of intersection (-118, -93, 44, 38)
NR.W.AQS.d.bound <- matrix(c(abs(38 - NR.W.AQS$Latitude), 
                             abs(44 - NR.W.AQS$Latitude)), ncol = 2)

NR.W.DS.d.bound <- matrix(c(abs(38 - NR.W.DS$Latitude), 
                            abs(44 - NR.W.DS$Latitude)), ncol = 2)

estimates.6 <- coef(mle(likelihoodFunRV2, start = list(a1 = 0, a2 = 0),
                        fixed = list(dist.cen = NR.W.AQS.d.cent, Y = NR.W.AQS$yi, 
                                     dist = NR.W.AQS.d.bound, mu = NR.W.AQS.mu, 
                                     sd = NR.W.AQS.sd)))

spliced.NR.W <- smoothEstimate2(a1 = estimates.6[1], a2 = estimates.6[2], 
                                dist.cen = NR.W.DS.d.cent, 
                                dist = NR.W.DS.d.bound, mu = NR.W.DS.mu)


###################
##Validation NW, NR

set.seed(1234)

test <- sample(45, 10)

estimates.3.cv <- coef(mle(likelihoodFunRV2, start = list(a1 = 0.1, a2 = 0.2), 
                        fixed = list(Y = NW.NR.AQS2[-test], 
                                     dist = NW.NR.AQS.d[-test, ], 
                                     mu = NW.NR.AQS.mu[-test, ], sd = NW.NR.AQS.sd[-test, ], 
                                     dist.cen = NW.NR.AQS.cent[-test])))

spliced.NW.NR.RV.new.cv <- smoothEstimate2(a1 = estimates.3.cv[1], a2 = estimates.3.cv[2],
                                        dist.cen = DS.d.cent.lon,
                                        dist = DS.d.bound, mu = DS.y)

NW.NR.matched <- match(AQS.grid$Loc_Label1, NW.NR.DS$Loc_Label1)

spliced.NW.NR.RV.new.cv <- (spliced.NW.NR.RV.new.cv[NW.NR.matched])[test]

#MSE from spliced data on test set
mean((NW.NR.AQS2[test] - spliced.NW.NR.RV.new.cv)^2)

#MSE from spliced data an 45 DS sites
spliced.NW.NR.RV.new.DS <- spliced.NW.NR.RV.new[NW.NR.matched]
mean((AQS.grid$yi - spliced.NW.NR.RV.new.DS)^2)

#MSE using all points on those selected by 'test'
mean((NW.NR.AQS2[test] - spliced.NW.NR.RV.new.DS[test])^2)

#MSE using NW DS estimates
mean((NW.NR.AQS2 - NW.NR.AQS.mu[, 1])^2)

#MSE using NR DS estimates
mean((NW.NR.AQS2 - NW.NR.AQS.mu[, 2])^2)

#Read National Estimates, keep predictions for AQS sites
national.DS <- read.csv("data/Results/pm25_2014_National/RESULTS/Predictions.csv")
national.DS <- national.DS[match(AQS.grid$Loc_Label1, national.DS$Loc_Label1), ]

#MSE for national vs AQS sites
mean((AQS.grid$yi - national.DS$Prediction)^2)

#Read IMPROVE readings, keep readings in the intersection
NW.NR.improve <- read.csv("data/NW_NR_improve.csv")
spliced.improve <- merge(NW.NR.improve, NW.NR.DS, by.x = "DS_Lat", by.y = "Latitude")

#MSE for spliced vs IMPROVE
mean((spliced.improve$spliced - spliced.improve$PM25_value)^2)

#MSE for DS NW vs IMPROVE
mean((spliced.improve$Prediction - spliced.improve$PM25_value)^2)

#MSE for DS NR vs IMPROVE
mean((spliced.improve$Prediction.1 - spliced.improve$PM25_value)^2)




setwd('C://Users//jeongwha//Desktop//EPA_SAMSI')



source("R/functions.R")

AQS.grid <- read.csv("data/AQS_DATA.csv")
NW.NR <- read.csv("data/combinedOfNW_NR.csv")

NW.NR.AQS2 <- AQS.grid$yi
NW.NR.AQS.mu <- matrix(c(NW.NR$m1, NW.NR$m2), ncol = 2)
NW.NR.AQS.sd <- matrix(c(NW.NR$s1, NW.NR$s2), ncol = 2)
NW.NR.AQS.d <- matrix(c(NW.NR$d.to.R_1.edge, NW.NR$d.to.R_2.edge), ncol = 2)

##likelihood function for new method : Y : vector(n by 1),  dist, mu, sd :matrix(n by 2) 

likelihoodFunction2<-function(phi, Y, dist, mu, sd) {     

 w <- exp(-phi*dist)
 new_mu<-apply(w*mu,1,sum)
  new_s<-apply((w^2)*(sd^2), 1,sum)
  new_s<-(new_s)^{1/2}
  likelihood <- prod(dnorm(Y, new_mu, new_s))

  return(-(log(likelihood))) }

  phi.mle2 <- (coef(mle(likelihoodFunction2, start = list(phi = .4), fixed = list(
  Y = NW.NR.AQS2, dist = NW.NR.AQS.d, mu = NW.NR.AQS.mu, sd = NW.NR.AQS.sd))))[1]
#phi.mle<-1.6
NW.NR.DS <- read.csv("data/combinedOverlap.csv")
#nrow(NW.NR.DS)
#DS.d <- matrix(c(118 + NW.NR.DS$Longitude, -1*NW.NR.DS$Longitude - 109), ncol = 2)
DS.y <- matrix(c(NW.NR.DS$Prediction, NW.NR.DS$Prediction.1), ncol = 2)
smoothEstimate1 <- function(d, mu,phi.mle2){
  #Returns combined estimates for each grid value in the intersection
  #
  #Args:
  #  d: Distance, measured in degrees of longitude, from each grid point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: All Downscaler estimates, Nx2 matrix.
  #
  #Returns:
  # Combined estimates for each grid point
  w <- exp(-1*phi.mle2*d)
  f <- (1/apply(w, 1, sum))*apply(w*mu, 1, sum)
}

spliced.NW.NR <- smoothEstimate1(DS.d, DS.y, phi.mle2)
#len<-length(spliced.NW.NR)
#NW.NR.DS$Latitude

#NW.NR.DS$Longitude

#howManyPredictions <- dim(Predictions)[1]

DATA.NR<-read.csv("data/Predictions.csv")
DATA.NR.Q1 <-DATA.NR[DATA.NR$Date=="Jan-01-2014",] 


  # take first prediction's location
  # find its location in DATA.Q1usa dataset
     newPredictionLocation <- match(NW.NR.DS$Loc_Label1,DATA.NR.Q1$Loc_Label1 )
  #replace the found location with new prediction
           DATA.NR.Q1$Prediction[newPredictionLocation]<- spliced.NW.NR


  setwd('C://Users//jeongwha//Desktop')
MapDATA <- DATA.NR.Q1
plotname <- "NationalPred.tiff"
plottitle <- "PM2.5 Seasonal Average: National Q1"
colorlabel <- "PM2.5 Seasonal Ave"



#mypalPred <- colorRampPalette(c("blue2",  "green4", "yellow1", "red4", "purple3")) 
#mypalPred attempts to scale colors to match Python output.  May need to play with
mypalPred <- colorRampPalette(c("blue2","blue1", "deepskyblue", "darkturquoise", "green3", "yellow1", "orange", "orangered", "red", "red2", "red3", "red4", "darkred")) #, "purple"
mypalStDev <- mypalPred
mypalDiff <- colorRampPalette(c("blue2",  "#ffe6e6", "red")) 

#plot all states with ggplot
#qplot(x=MapDATA$Lon, y=MapDATA$Lat, color=MapDATA$Prediction)
all_states <- map_data("state") #creates state boundaries
Sys.time()
p <- ggplot()
p <- p + geom_polygon( data=all_states, aes(x=long, y=lat, group = group),colour="white", fill="grey50" ) #grey10 is almost black
#Add DS output, small pointsize otherwise overlap
p <- p + geom_point(data=MapDATA, aes(x=MapDATA$Lon, y=MapDATA$Lat, color=MapDATA$Prediction), size=.075)  +
  scale_colour_gradientn(colours = mypalPred(8)) + theme_bw()
#
p <- p +  coord_map("cylindrical") + labs(colour = colorlabel) 
#Add AQS sites
p <- p + geom_point(data=DATA.aqs, aes(x=DATA.aqs$Lon, y=DATA.aqs$Lat), color="black", size=.1)
#overlay state boundaries; no fill
p <- p + geom_polygon( data=all_states, aes(x=long, y=lat, group = group),colour="white", fill="NA",size=.2 ) 

#p <- p + ggtitle(plottitle)
###DON'T OUTPUT THIS GRAPH IN RSTUDIO!
######p #takes hour(s) to output when mapping all DS Predictions

Sys.time() #5-10 sec
res.pic<-600
ggsave(plotname,dpi=res.pic)

Sys.time()         
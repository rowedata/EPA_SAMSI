#CREATE MAPS FOR NW NR REGIONS

library(rgdal) 
library(sp)
library(ggplot2)
library(geosphere)
library(ggmap)
library(data.table)

data.dir <- "D:/Ohio State/SAMSI IMSM/EPA/Data"  
cmaq.dir <- "D:/Ohio State/SAMSI IMSM/EPA/Data/DSinput/ForAdam"
aqs.dir <- "D:/Ohio State/SAMSI IMSM/EPA/Data/DSinput/ForAdam"  
ds.dir <- "D:/Ohio State/SAMSI IMSM/EPA/Data/Results"
export.lib <- "D:/Ohio State/SAMSI IMSM/EPA_SAMSI/figs"
res.pic <- 600

Cb.folder <- "/cb_2017_us_nation_5m"
Cb.shpfile <- "cb_2017_us_nation_5m"
regional.dir <- "pm25_2014_"
regions <- c("NE", "NR", "NW", "OV", "S", "SE", "SW", "UM", "W")
ds.output <- "RESULTS/Predictions.csv"


#Reading data

setwd(aqs.dir) 
DATA.aqs <- read.csv("ds.input.aqs.pm25.2014.quarterly.csv", header=TRUE)
head(DATA.aqs) #4410 obs of 6 vars
#Site POC       Date      Lat        Lon      Conc
DATA.aqs <- data.frame(DATA.aqs)
setwd(ds.dir) 

DATA <- read.csv(paste0(regional.dir,"National/", ds.output), header=TRUE) # NR only


# LOAD REGIONAL DOWNSCLER PREDICTIONS 
DATA.NR <- read.csv(paste0(regional.dir,"NR/", ds.output), header=TRUE) # NR only
DATA.NW <- read.csv(paste0(regional.dir,"NW/", ds.output), header=TRUE) # NW only
DATA.W <- read.csv(paste0(regional.dir,"W/", ds.output), header=TRUE) # W only

DATA.NR <- data.frame(DATA.NR)
DATA.NW <- data.frame(DATA.NW)
DATA.W <- data.frame(DATA.W)

DATA.NR.Q1 <-DATA.NR[DATA.NR$Date=="Jan-01-2014",] 
DATA.NW.Q1 <-DATA.NW[DATA.NW$Date=="Jan-01-2014",] 
DATA.W.Q1 <-DATA.W[DATA.W$Date=="Jan-01-2014",] 

#Date Loc_Label1 Latitude Longitude Prediction SEpred
DATA.Q1 <- DATA[DATA$Date=="Jan-01-2014",] #137241
DATA.Q2 <- DATA[DATA$Date=="Jan-02-2014",] #137241
DATA.Q3 <- DATA[DATA$Date=="Jan-03-2014",] #137241
DATA.Q4 <- DATA[DATA$Date=="Jan-04-2014",] #137241

#Use Censur Bureau shapefile
#https://www.census.gov/geo/maps-data/data/cbf/cbf_nation.html 
us_cb.shp<-readOGR(paste0(data.dir,Cb.folder), Cb.shpfile) 
us_cb.shp<-spTransform(us_cb.shp,CRS("+proj=longlat"))
class(us_cb.shp)
####

#Called "Nat" because this is the National file (which includes some ocean/Canada)
# Nat.pts_proj <- SpatialPoints(cbind(DATA.Q1$Longitude,DATA.Q1$Latitude),
#                               proj4string=CRS("+proj=longlat"))
# 
# keep_cb_proj <- over(Nat.pts_proj, us_cb.shp)
# DATA.Q1["InConus"]<- is.finite(keep_cb_proj$GEOID)
# table(DATA.Q1$InConus)
# DATA.Q1usa <- DATA.Q1[DATA.Q1$InConus,]

#############################################################################
setwd("D:/Ohio State/SAMSI IMSM/EPA_SAMSI/data")
  
# TRY TO REPLACE NR PREDICTIONS WITH OVERLAP NR
head(DATA.NW.Q1)
dim(DATA.NW.Q1)

predictions <- read.csv("Overlap.csv", header=TRUE)
howManyPredictions <- dim(predictions)[1]

for (i in 1:howManyPredictions){
  # find its location in DATA.Q1usa dataset
  newPredictionLocation <- match(predictions$Loc_Label1[i], DATA.NW.Q1$Loc_Label1)
  #replace the found location with new prediction
  DATA.NW.Q1$Prediction[newPredictionLocation] = spliced.NW.W[i]
}

NW.matched <- match(predictions$Loc_Label1, DATA.W.Q1$Loc_Label1)
DATA.NW.Q1$Prediction[NW.matched] = spliced.NW.W

######################################################################################
# REPLACE NR OVERLAP with overlap predictions

head(DATA.NW.Q1)
dim(DATA.NW.Q1)
howManyPredictions <- dim(predictions)[1]

for (i in 1: howManyPredictions) {
  predictions$Loc_Label1[i] # take first prediction's location
  # find its location in DATA.Q1usa dataset
  newPredictionLocation <- match(predictions$Loc_Label1[i],DATA.NR.Q1$Loc_Label1 )
  #replace the found location with new prediction
  DATA.NR.Q1$Prediction[newPredictionLocation] = spliced.NW.NR.RV.4[i]
}


setwd(export.lib)

MapDATA <- DATA.NR.Q1
plotname <- "NationalPredModel4new.tiff"
plottitle <- "PM2.5 Seasonal Average: National Q1"
colorlabel <- "PM2.5 Seasonal Ave"

#mypalPred <- colorRampPalette(c("blue2",  "green4", "yellow1", "red4", "purple3")) 
#mypalPred attempts to scale colors to match Python output.  May need to play with
mypalPred <- colorRampPalette(c("blue2","blue1", "deepskyblue", "darkturquoise", "green3", "yellow1", "orange", "orangered", "red", "red2", "red3", "red4", "darkred")) #, "purple"
mypalStDev <- mypalPred
mypalDiff <- colorRampPalette(c("blue2",  "#ffe6e6", "red")) 

#plot all states with ggplot
all_states <- map_data("state") #creates state boundaries
Sys.time()
  p <- ggplot()
  p <- p + geom_polygon( data=all_states, aes(x=long, y=lat, group = group),colour="white", fill="grey50" ) #grey10 is almost black

#Add DS output, small pointsize otherwise overlap
  p <- p + geom_point(data=MapDATA, aes(x=MapDATA$Lon, y=MapDATA$Lat, color=MapDATA$Prediction), size=.075)  +
    scale_colour_gradientn(colours = mypalPred(8)) + theme_bw()
  p <- p +  coord_map("cylindrical") + labs(colour = colorlabel) 

#Add AQS sites
  p <- p + geom_point(data=DATA.aqs, aes(x=DATA.aqs$Lon, y=DATA.aqs$Lat), color="black", size=.1)
#overlay state boundaries; no fill
  p <- p + geom_polygon( data=all_states, aes(x=long, y=lat, group = group),colour="white", fill="NA",size=.2 ) 

#p <- p + ggtitle(plottitle)

ggsave(plotname,dpi=res.pic)








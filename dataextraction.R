rm(list = ls())
library(readr)






###########################
# function: extracts data in the overlapping zone, find AQS sites #
###########################
AQSsites <- function(long_l, long_r, lat_d, lat_u, datavals){
  datavals[which(datavals$Lon > long_l & datavals$Lon < long_r & datavals$Lat > lat_d & datavals$Lat < lat_u), ] 
}

DSsites <- function(long_l, long_r, lat_d, lat_u, datavals){
  datavals[which(datavals$Longitude > long_l & datavals$Longitude < long_r & datavals$Latitude > lat_d & datavals$Latitude < lat_u), ] 
}

#####AQS_overlap_sites <- sites(-118, -109, 40, 51, NW_AQS)


###########################
# function: extracts data in the overlapping zone by quarter #
###########################
AQSquart <- function(datavals, quarter){
  if (quarter == 1){
    q <- "2014-01-01"
  }else if (quarter ==2){
    q <- "2014-01-02"
  }else if (quarter ==3){
    q <- "2014-01-03"
  }else if (quarter ==4){
    q <- "2014-01-04"
  }else {
    print("ERROR")
  }
  datavals[datavals$Date == q, ]
}

DSquart <- function(datavals, quarter){
  if (quarter == 1){
    q <- "Jan-01-2014"
  }else if (quarter ==2){
    q <- "Jan-02-2014"
  }else if (quarter ==3){
    q <- "Jan-03-2014"
  }else if (quarter ==4){
    q <- "Jan-04-2014"
  }else {
    print("ERROR")
  }
  datavals[datavals$Date == q, ]
}
#breaks up AQS data for each quarter
######AQS_overlap_sites.Q1 <- quart(AQS_overlap_sites,1) #45
######AQS_overlap_sites.Q2 <- quart(AQS_overlap_sites,2) #37
######AQS_overlap_sites.Q3 <- quart(AQS_overlap_sites,3) #42
######AQS_overlap_sites.Q4 <- quart(AQS_overlap_sites,4) #33


###########################
# function: calculates the vector 2-norm #
###########################
norm_vec <- function(x) sqrt(sum(x^2))


###########################
# function: calculates the DS data for AQS sites in overlap area  #
###########################
grid <- function(AQS,DS_R1, DS_R2){
  ns <- nrow(AQS) #45
  nds <- nrow(DS_R1) #5978
  
  #sets inital values for distance and minimum row
  d_old <- 1000000
  min_grid <- 0
  
  #takes mu and SE from the NR "zone 2" region
  zone2 <- cbind(DS_R2$Prediction, DS_R2$SEpred)
  name <- c("mu_2", "SE_2")
  colnames(zone2) <- name
  
  #initializes AQS_grid
  grid <- cbind(DS_R1[1:ns, ], zone2[1:ns, ])
  colnames(grid)[colnames(grid) == "Prediction"] <- "mu_1"
  colnames(grid)[colnames(grid) == "SEpred"] <- "SE_1"
  
  #calculates where the AQS locations with respect to the grid points
  for (j in 1:ns){
    for (i in 1:nds){
      
      #for a fixed AQS location, it calculates the norm of the difference to all of the grid points in the overlap region
      vec <- c(DS_R1$Longitude[i] - AQS$Lon[j], DS_R1$Latitude[i] - AQS$Lat[j])
      d_new <- norm_vec(vec) 
      
      #keeps track of the minimum distance and corresponding row value
      if (d_new<d_old){
        d_old <- d_new
        min_grid <- i
      }
    }
    #print(min_grid)
    
    #assigns the NW and NR DS data for the gridpoint the AQS "station" is located to AQS_grid
    grid[j, ] <- cbind(DS_R1[min_grid, ], zone2[min_grid,1], zone2[min_grid,2])
    
    #reassigns initial values for the distance and minimum row
    d_old <- 1000000
    min_grid <- 0
  }
  
  return(grid)
  
}

#####AQS_grid <- grid(AQS_overlap_sites.Q1, DS_overlap_sites_NW.Q1, DS_overlap_sites_NR.Q1)
  

###########################
# function: calculates the distance to the zones #
###########################
distance <- function(latvslong, AQS, long_l, long_r, lat_d, lat_u, AQS_grid){
  #calulate distance to edge of boundaries
  
  ns <- nrow(AQS) #45
  
  #initializes distance matrix
  distance <- matrix(0, ns, 2)
  
  if (latvslong == "lat"){
    #height of the overlap area
    a <- abs(lat_d - lat_u)
    
    #calulates distance to each boundary
    for (i in 1:ns){
      #takes row in ith AQS_grid
      rowtable <- AQS_grid[i, ]
      #calculates distance to bottom boundary (edge of region 1)
      d1 <- abs(rowtable$Lattitude - (lat_d))
      #calculates distance to top boundary (edge of region 2)
      d2 <- abs(rowtable$Lattitude - (lat_u))
      #calulates the difference between the width of the overlapping region and the distance to the boundary
      #we needed to do this to ensure that the further into a region a point is, the larger the weight of that mean is
      distance[i, ] <- c(a - d1, a - d2)
    }
  } else if (latvslong == "long"){
    #width of the overlap area
    a <- abs(long_r - long_l) #9
    
    #calulates distance to each boundary
    for (i in 1:ns){
      #takes row in ith AQS_grid
      rowtable <- AQS_grid[i, ]
      #calculates distance to left boundary (edge of region 1)
      d1 <- abs(rowtable$Longitude - (long_r))
      #calculates distance to right boundary (edge of region 2)
      d2 <- abs(rowtable$Longitude - (long_l))
      #calulates the difference between the width of the overlapping region and the distance to the boundary
      #we needed to do this to ensure that the further into a region a point is, the larger the weight of that mean is
      distance[i, ] <- c(a - d1, a - d2)
    }
  }
  #adds the distance to region 1 and region 2 to the table AQS_grid
  name <- c("d to R_1 edge", "d to R_2 edge")
  colnames(distance) <- name
  AQS_grid <- cbind(AQS_grid, distance)
  
  #adds AQS data to the table AQS_grid
  yi <- AQS$Conc
  AQS_grid <- cbind(AQS_grid,yi)
  
}

#####AQS_grid2 <- distance("long", AQS_overlap_sites.Q1, -118, -109, 40, 51, AQS_grid)
  


NW_AQS <- read_csv("~/Documents/IMSM/Working progress/NW_AQS.csv")
NW_DS <- read_csv("~/Documents/IMSM/Working progress/NW_DS.csv")
NR_DS <- read_csv("~/Documents/IMSM/Working progress/NR_DS.csv")

organized_data <- function(long_l, long_r, lat_d, lat_u, AQS, DS_R1, DS_R2, latvslong){
  #AQS = NW_AQS
  AQS_overlap_sites <- AQSsites(long_l, long_r, lat_d, lat_u, AQS)
  AQS_overlap_sites.Q1 <- AQSquart(AQS_overlap_sites,1) #45
  AQS_overlap_sites.Q2 <- AQSquart(AQS_overlap_sites,2) #37
  AQS_overlap_sites.Q3 <- AQSquart(AQS_overlap_sites,3) #42
  AQS_overlap_sites.Q4 <- AQSquart(AQS_overlap_sites,4) #33
  
  #DS_R1 = NW_DS
  DS_R1_overlap_sites <- DSsites(long_l, long_r, lat_d, lat_u, DS_R1)
  DS_R1_overlap_sites.Q1 <- DSquart(DS_R1_overlap_sites,1) #45
  DS_R1_overlap_sites.Q2 <- DSquart(DS_R1_overlap_sites,2) #37
  DS_R1_overlap_sites.Q3 <- DSquart(DS_R1_overlap_sites,3) #42
  DS_R1_overlap_sites.Q4 <- DSquart(DS_R1_overlap_sites,4) #33
  
  #DS_R2 = NR_DS
  DS_R2_overlap_sites <- DSsites(long_l, long_r, lat_d, lat_u, DS_R2)
  DS_R2_overlap_sites.Q1 <- DSquart(DS_R2_overlap_sites,1) #45
  DS_R2_overlap_sites.Q2 <- DSquart(DS_R2_overlap_sites,2) #37
  DS_R2_overlap_sites.Q3 <- DSquart(DS_R2_overlap_sites,3) #42
  DS_R2_overlap_sites.Q4 <- DSquart(DS_R2_overlap_sites,4) #33
  
  AQS_grid <- grid(AQS_overlap_sites.Q1, DS_R1_overlap_sites.Q1, DS_R2_overlap_sites.Q1)
  AQS_grid2 <- distance("long", AQS_overlap_sites.Q1, long_l, long_r, lat_d, lat_u, AQS_grid)
  return(AQS_grid2)
}

organized_data(-118, -109, 40, 51, NW_AQS, NW_DS, NR_DS, "long")







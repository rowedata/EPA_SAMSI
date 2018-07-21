###########################
# function: converts region boundaries to overlap boundaries #
#longlat must be formatted as (long_l, long_r, lat_u, lat_d) 
###########################
bound <- function(ll_R1, ll_R2){
  b_ll <- max(ll_R1[1], ll_R2[1])
  b_lr <- min(ll_R1[2], ll_R2[2])
  b_ld <- max(ll_R1[4], ll_R2[4])
  b_lu <- min(ll_R1[3], ll_R2[3])
  
  if(b_lr <= b_ll || b_ld >= b_lu){
    return("ERROR: NO OVERLAP!!!!!!")
  }else{
    return(c(b_ll, b_lr, b_ld, b_lu))
  }
  
}

###########################
# function: extracts DS data in the overlapping zone #
###########################
DSsites <- function(long_l, long_r, lat_d, lat_u, datavals){
  datavals[which(datavals$Longitude > long_l & datavals$Longitude < long_r & datavals$Latitude > lat_d & datavals$Latitude < lat_u), ] 
}

###########################
# function: breaks up DS data according to quarter #
###########################
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


###########################
# function: gets  DS data for datavals for a specific 
  #quarter in the overlap of regions governed by ll_R1, ll_R2 #
###########################
getdata <- function(ll_R1, ll_R2, datavals, quarter){
  #input: ll_R1 <- (long_l, long_r, lat_u, lat_d) of region 1
         #ll_R2 <- (long_l, long_r, lat_u, lat_d)  of region 2
         #datavals <- region you want the values of (i.e. NW_DS)
         #quarter <- quarter you want (i.e. 1, 2, 3, or 4)
  #output: DS data for quarter in overlap region of region 1 and region 2
  bnd <- bound(ll_R1, ll_R2)
  
  long_l <- bnd[1]
  long_r <- bnd[2]
  lat_u <- bnd[3]
  lat_d <- bnd[4]
  
  sites <- DSsites(long_l, long_r, lat_d, lat_u, datavals)
  sites.Q1 <- DSquart(datavals, quarter)
  
  return(sites.Q1)
  
}
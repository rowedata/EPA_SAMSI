
start.time = Sys.time()

###############################################################################################

#################
### FUNCTIONS ###
#################

#=== Function 1: Extracts data in the overlapping zone for AQS ===#
AQSsites = function(long_l, long_r, lat_d, lat_u, datavals){
  
  subset(datavals, (Lat < lat_u & Lat > lat_d & Lon > long_l  & Lon < long_r))
  
}


#=== Function 2: Extracts data in the overlapping zone for DS ==#
DSsites = function(long_l, long_r, lat_d, lat_u, datavals){
 
  subset(datavals, (Latitude < lat_u & Latitude > lat_d & Longitude > long_l  & Longitude < long_r))
  
}


#=== Function 3: extracts data in the overlapping zone by quarter (AQS)#
AQSquart = function(datavals, quarter){
  if (quarter == 1){
    
    q = "2014-01-01"
    
  }else if (quarter == 2){
    
    q = "2014-01-02"
    
  }else if (quarter == 3){
    
    q = "2014-01-03"
  }else if (quarter == 4){
    
    q = "2014-01-04"
  }else {
    print("ERROR")
  }
  datavals[datavals$Date == q, ]
}


#=== Function 4: extracts data in the overlapping zone by quarter (DS)#
DSquart = function(datavals, quarter){
  
  if (quarter == 1){
    
    q = "Jan-01-2014"
  }else if (quarter == 2){
    
    q = "Jan-02-2014"
  }else if (quarter == 3){
    
    q = "Jan-03-2014"
  }else if (quarter == 4){
    
    q = "Jan-04-2014"
  }else {
    print("ERROR")
  }
  
  datavals[datavals$Date == q, ]
  
}


#=== Function 5: Calculates the vector 2-norm ===#
norm_vec = function(x) sqrt(sum(x^2))


#=== Function 6: Calculates the DS data for AQS sites in overlap area  ===#
grid = function(AQS, DS_R1, DS_R2){
  
  ns  = nrow(AQS) #45
  nds = nrow(DS_R1) #5978
  
  #=== sets inital values for distance and minimum row ===#
  d_old = 1000000
  min_grid = 0
  
  # takes mu and SE from the NR "zone 2" region #
  zone2 = cbind(DS_R2$Prediction, DS_R2$SEpred)
  colnames(zone2) = c("mu_2", "SE_2")
  
  # initializes AQS_grid #
  grid = cbind(DS_R1[1:ns, ], zone2[1:ns, ])
  colnames(grid)[colnames(grid) == "Prediction"] = "mu_1"
  colnames(grid)[colnames(grid) == "SEpred"] = "SE_1"
  
  # calculates where the AQS locations with respect to the grid points #
  for (j in 1:ns){
    for (i in 1:nds){
      
      # for a fixed AQS location, it calculates the norm of the difference to all of the grid points in the overlap region
      vec   = c(DS_R1$Longitude[i] - AQS$Lon[j], DS_R1$Latitude[i] - AQS$Lat[j])
      d_new = norm_vec(vec) 
      
      #keeps track of the minimum distance and corresponding row value
      if (d_new < d_old){
          d_old  = d_new
        min_grid = i
      }
    }
    #print(min_grid)
    
    #assigns the NW and NR DS data for the gridpoint the AQS "station" is located to AQS_grid
    grid[j, ] = cbind(DS_R1[min_grid, ], zone2[min_grid, 1], zone2[min_grid, 2])
    
    #reassigns initial values for the distance and minimum row
    d_old    = 1000000
    min_grid = 0
  }
  
  return(grid)
  
}


#=== Function 7: calculates the distance to the zones ===#
distance = function(latvslong, AQS, long_l, long_r, lat_d, lat_u, AQS_grid){
 
   # calulate distance to edge of boundaries #
  ns = nrow(AQS) #45
  
  #initializes distance matrix
  distance = matrix(0, ns, 2)
  
  if (latvslong == "lat"){
    
    # height of the overlap area #
    a = abs(lat_d - lat_u)
    
    #calulates distance to each boundary
    for (i in 1:ns){
      
      rowtable = AQS_grid[i, ]
      
      # calculates distance to bottom boundary (edge of region 1) #
      d1 = abs(rowtable$Lattitude - (lat_d))
      d2 = abs(rowtable$Lattitude - (lat_u))
      
      distance[i, ] = c(a - d1, a - d2)
    }
    
  } else if (latvslong == "long"){
    
    # width of the overlap area #
    a = abs(long_r - long_l) 
    
    # calulates distance to each boundary #
    for (i in 1:ns){
      
      rowtable = AQS_grid[i, ]
      
      d1 = abs(rowtable$Longitude - (long_r))
      d2 = abs(rowtable$Longitude - (long_l))
      
      distance[i, ] = c(a - d1, a - d2)
      
    }
  }
  # adds the distance to region 1 and region 2 to the table AQS_grid #
  name = c("d1", "d2")
  colnames(distance) = name
  AQS_grid = cbind(AQS_grid, distance)
  
  Y_i = AQS$Conc
  AQS_grid = cbind(AQS_grid, Y_i)
  
}

#=== Function 8: Organized data ===#
organized_data = function(ll_R_1, ll_R_2, AQS, DS_R1, DS_R2, latvslong){
  
  bnd = bound(ll_R_1, ll_R_2)
  
  if(bnd[[1]] == "ERROR: NO OVERLAP!!!!!!"){
    
    return("ERROR: NO OVERLAP!!!!!!")
    
  }else{
    long_l = bnd[[1]]
    long_r = bnd[[2]]
    lat_d  = bnd[[3]]
    lat_u  = bnd[[4]]
  
  # AQS = NW_AQS #
  AQS_overlap_sites    = AQSsites(long_l, long_r, lat_d, lat_u, AQS)
  AQS_overlap_sites.Q1 = AQSquart(AQS_overlap_sites, 1) #45
  AQS_overlap_sites.Q2 = AQSquart(AQS_overlap_sites, 2) #37
  AQS_overlap_sites.Q3 = AQSquart(AQS_overlap_sites, 3) #42
  AQS_overlap_sites.Q4 = AQSquart(AQS_overlap_sites, 4) #33
  
  # DS_R1 = NW_DS #
  DS_R1_overlap_sites    = DSsites(long_l, long_r, lat_d, lat_u, DS_R1)
  DS_R1_overlap_sites.Q1 = DSquart(DS_R1_overlap_sites, 1) #45
  DS_R1_overlap_sites.Q2 = DSquart(DS_R1_overlap_sites, 2) #37
  DS_R1_overlap_sites.Q3 = DSquart(DS_R1_overlap_sites, 3) #42
  DS_R1_overlap_sites.Q4 = DSquart(DS_R1_overlap_sites, 4) #33
  
  # DS_R2 = NR_DS #
  DS_R2_overlap_sites    = DSsites(long_l, long_r, lat_d, lat_u, DS_R2)
  DS_R2_overlap_sites.Q1 = DSquart(DS_R2_overlap_sites, 1) #45
  DS_R2_overlap_sites.Q2 = DSquart(DS_R2_overlap_sites, 2) #37
  DS_R2_overlap_sites.Q3 = DSquart(DS_R2_overlap_sites, 3) #42
  DS_R2_overlap_sites.Q4 = DSquart(DS_R2_overlap_sites, 4) #33
  
  AQS_grid  = grid(AQS_overlap_sites.Q1, DS_R1_overlap_sites.Q1, DS_R2_overlap_sites.Q1)
  AQS_grid2 = distance("long", AQS_overlap_sites.Q1, long_l, long_r, lat_d, lat_u, AQS_grid)
  return(AQS_grid2)
  
  }
  
}  


#=== Function 9: Bound ===#
bound = function(ll_R1, ll_R2){
  
  b_ll = max(ll_R1[1], ll_R2[1])
  b_lr = min(ll_R1[2], ll_R2[2])
  b_ld = max(ll_R1[4], ll_R2[4])
  b_lu = min(ll_R1[3], ll_R2[3])
  
  if(b_lr <= b_ll || b_ld >= b_lu){
    return("ERROR: NO OVERLAP!!!!!!")
    
  }else{
    return(c(b_ll, b_lr, b_ld, b_lu))
  }
  
}

#=== Function 10: Optimzation function ===#

Likelihood_func  =  function(phi){
  
  w =  exp(-phi*d)
  
  R = prod((1/(apply(w, 1, sum)))*(apply(w*dnorm(Y, mu, sd), 1, sum)))
  
  -(log(R))
}


###############################################################################################

#############################
#=== WORKING DIRECTORIES ===#
#############################

nw_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_NW/RESULTS"
w_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_W/RESULTS"
nr_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_NR/RESULTS"
sw_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_SW/RESULTS"
um_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_UM/RESULTS"
s_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_S/RESULTS"
ov_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_OV/RESULTS"
ne_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_NE/RESULTS"
se_downscaler.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/Results/pm25_2014_SE/RESULTS"


aqs.dir = "C:/Users/dofori/Desktop/SAMSI/EPA/Data/DSinput/ForAdam"
Main.dir = "C:/Users/dofori/Desktop/SAMSI/Results_folder"

###############################################################################################

##################
#=== AQS DATA ===#
##################

setwd(aqs.dir)
NW_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.NW.csv")

setwd(aqs.dir)
W_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.W.csv")

setwd(aqs.dir)
NR_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.NR.csv")

setwd(aqs.dir)
SW_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.SW.csv")

setwd(aqs.dir)
UM_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.UM.csv")

setwd(aqs.dir)
S_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.S.csv")

setwd(aqs.dir)
OV_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.OV.csv")

setwd(aqs.dir)
NE_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.NE.csv")

setwd(aqs.dir)
SE_AQS = read.csv("ds.input.aqs.pm25.2014.quarterly.SE.csv")

###############################################################################################

#################
#=== DS DATA ===#
#################

setwd(nw_downscaler.dir)
NW_DS = read.csv("Predictions.csv")

setwd(w_downscaler.dir)
W_DS = read.csv("Predictions.csv")

setwd(nr_downscaler.dir)
NR_DS = read.csv("Predictions.csv")

setwd(sw_downscaler.dir)
SW_DS = read.csv("Predictions.csv")

setwd(um_downscaler.dir)
UM_DS = read.csv("Predictions.csv")

setwd(s_downscaler.dir)
S_DS = read.csv("Predictions.csv")

setwd(ov_downscaler.dir)
OV_DS = read.csv("Predictions.csv")

setwd(ne_downscaler.dir)
NE_DS = read.csv("Predictions.csv")

setwd(se_downscaler.dir)
SE_DS = read.csv("Predictions.csv")

###############################################################################################

####################################
#==== DATASET EXTRACTION (SOME)====#
####################################

#== S_SE ==#
latlong1 = c(-108, -86, 42, 24)
latlong2 = c(-88, -73, 42, 23)

organised.list = organized_data(latlong1, latlong2, S_AQS, S_DS, SE_DS, "long")

#== NW_NR ==#
latlong1 = c(-127, -109, 51, 40)
latlong2 = c(-118, -93, 51, 38)

organised.list2 = organized_data(latlong1, latlong2, NW_AQS, NW_DS, NR_DS, "long")

#== OV_SE ==#
latlong1 = c(-97, -76, 45, 33)
latlong2 = c(-88, -73, 42, 23)

organised.list3 = organized_data(latlong1, latlong2, OV_AQS, OV_DS, SE_DS, "long")

#==SW_S ==#
latlong1 = c(-116, -100, 44, 29)
latlong2 = c(-108, -86, 42, 24)

organised.list4 = organized_data(latlong1, latlong2, SW_AQS, SW_DS, S_DS, "long")


latlong1 = c(-99, -81, 51, 39)
latlong2 = c(-97, -76, 45, 33)

organised.list5 = organized_data(latlong1, latlong2, UM_AQS, UM_DS, OV_DS, "lat")

##############################################################################################

################################
#=== EXPORT RESULTS TO FILE ===#
################################

setwd(Main.dir)
write.csv(organised.list, "S_n_SE.csv")
write.csv(organised.list2, "NW_n_NR.csv")
write.csv(organised.list3, "OV_n_SE.csv")
write.csv(organised.list4, "SW_n_S.csv")
write.csv(organised.list5, "UM_n_OV.csv")


#############################################################################################

######################################################
#=== Optimisation of phi value in Model averaging ===#
######################################################
 
t = 0

data1 = organised.list

t = t + 1
Y   =  matrix(c(data1$Y_i , data1$Y_i), ncol = 2)
mu  =  matrix(c(data1$mu_1, data1$mu_2), ncol = 2)
sd  =  matrix(c(data1$SE_1, data1$SE_1), ncol = 2)
d   =  matrix(c(data1$d1  , data1$d2  ), ncol = 2)


phi_optim[t] = (coef(mle(Likelihood_func, start = list(phi = 1))))[[1]]

###############################################################################################

data2 = organised.list2

t = t + 1
Y   =  matrix(c(data1$Y_i , data1$Y_i), ncol = 2)
mu  =  matrix(c(data1$mu_1, data1$mu_2), ncol = 2)
sd  =  matrix(c(data1$SE_1, data1$SE_1), ncol = 2)
d   =  matrix(c(data1$d1  , data1$d2  ), ncol = 2)

phi_optim[t] = (coef(mle(Likelihood_func, start = list(phi = 1))))[[1]]

###############################################################################################

data3 = organised.list3

t = t + 1
Y   =  matrix(c(data1$Y_i , data1$Y_i) , ncol = 2)
mu  =  matrix(c(data1$mu_1, data1$mu_2), ncol = 2)
sd  =  matrix(c(data1$SE_1, data1$SE_1), ncol = 2)
d   =  matrix(c(data1$d1  , data1$d2)  , ncol = 2)

phi_optim[t] = (coef(mle(Likelihood_func, start = list(phi = 1))))[[1]]

###############################################################################################

end.time = Sys.time()

time.taken = end.time - start.time



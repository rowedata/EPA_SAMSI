### include the AQS station data from the intersection of NR and NW
setwd('C://Users//sjhuv//Desktop//Misc//IMSM workshop//workshop material')
AQS_grid<-read.csv("AQS_DATA.csv")
Overlap<-read.csv('combinedOverlap.csv')
head(Overlap)
head(AQS_grid)
### the MLE of the estimated parameter
phi.NR.NW<-1.449767
weight<-function(distance,phi){
  return(exp(-phi*distance))
}

### bound of longitude of the overlap area
Overlap.long.upper.bound<-(-109)
Overlap.long.lower.bound<-(-118)

### for each grid cell, calculate the distance from it to two overlapping regions
grid.num<-nrow(Overlap)
weight.NW<-matrix(nrow = 1,ncol = grid.num)
weight.NR<-matrix(nrow = 1,ncol = grid.num)
### Calculate the weight based on the relative distances
for (i in 1:grid.num){
  ### define the two distances
  dist.NR<-abs(Overlap$Longitude[i]-Overlap.long.upper.bound)
  dist.NW<-abs(Overlap$Longitude[i]-Overlap.long.lower.bound)
  weight.NR[i]<-weight(dist.NR,phi.NR.NW)
  weight.NW[i]<-weight(dist.NW,phi.NR.NW)
}
weight.mean<-matrix(nrow = 1,ncol = grid.num)
weight.std<-matrix(nrow = 1,ncol = grid.num)
for (i in 1:grid.num){
  weight.sum<-weight.NR[i]+weight.NW[i]
  ### get the weighted mean
  weight.mean[i]<-weight.NR[i]/(weight.sum)*Overlap$Prediction.1[i]+
    weight.NW[i]/(weight.sum)*Overlap$Prediction[i]
  weight.std[i]<-sqrt((weight.NR[i]/weight.sum)*((Overlap$Prediction.1)^2+(Overlap$SEpred.1)^2)+
    (weight.NW[i]/weight.sum)*((Overlap$Prediction)^2+(Overlap$SEpred)^2)-(weight.mean[i])^2)
}
Overlap<-cbind(Overlap,t(weight.mean))
head(Overlap)
colnames(Overlap)[9]<-c('Weighted_Mean')
write.csv(Overlap,'Overlap.csv')

### Exploratory Data Analysis
### See the discrepancies between DS and Our Prediction
### Define the discrepancies
Relative.Discrepancy<-function(first.num,second.num){
  return(2*abs(first.num-second.num)/(first.num+second.num))
}
Discrepancy.NW<-matrix(nrow = 1,ncol = grid.num)
Discrepancy.NR<-matrix(nrow = 1,ncol = grid.num)
for (i in 1:grid.num){
  Discrepancy.NR[i]<-Relative.Discrepancy(Overlap$Prediction.1[i],Overlap$Weighted_Mean[i])
  Discrepancy.NW[i]<-Relative.Discrepancy(Overlap$Prediction[i],Overlap$Weighted_Mean[i])
}
### See the discrepancies between AQS readings and weighted mean in its grid
head(AQS_grid)
AQS.num<-nrow(AQS_grid)

### AQS_grid records the data for each AQS in intersection and its corresponding grid cell information
AQS.NW.weight<-matrix(nrow = 1,ncol = AQS.num)
AQS.NR.weight<-matrix(nrow = 1,ncol = AQS.num)
AQS.weighted.mean<-matrix(nrow = 1,ncol = AQS.num)
AQS.Discrepancy<-matrix(nrow = 1,ncol = AQS.num)
for (i in 1:AQS.num){
  AQS.NW.weight[i]<-weight(AQS_grid$d.to.R_1.edge[i],phi.NR.NW)
  AQS.NR.weight[i]<-weight(AQS_grid$d.to.R_2.edge[i],phi.NR.NW)
  weight.sum<-AQS.NW.weight[i]+AQS.NR.weight[i]
  AQS.weighted.mean[i]<-(AQS.NW.weight[i]/weight.sum)*(AQS_grid$mu_1[i])+AQS.NR.weight[i]/weight.sum*(AQS_grid$mu_2[i])
  AQS.Discrepancy[i]<-Relative.Discrepancy(AQS.weighted.mean[i],AQS_grid$yi[i])
}
hist(Discrepancy.NR)
hist(Discrepancy.NW)
hist(AQS.Discrepancy)
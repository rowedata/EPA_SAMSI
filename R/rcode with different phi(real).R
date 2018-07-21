
d<-seq(0,9,0.01);
gf<-function(phi)
{
exp(-phi*d)/(exp(-phi*d)+exp(-phi*(9-d)))
}
plot(d, gf(3), type="l", col="blue", xlab="distance",ylab="weight", 
     main="weight fucntion for different phi")

lines(d,gf(0), type="l", col="black")
lines(d, gf(1), type="l", col="red")

legend("topright", bty='n', xpd=NA,
       c("phi=0", "phi=1", "phi=3"), col=c("black", "red", "blue"), lty=c(1, 1))


#setwd("/R codes model/Model 3/0.5")
#source('variance.data.R')

proteins1<-data.matrix(read.csv("noise data 0.5%_r1.csv", row.names=1))
proteins2<-data.matrix(read.csv("noise data 0.5%_r2.csv", row.names=1))
proteins3<-data.matrix(read.csv("noise data 0.5%_r3.csv", row.names=1))
## Due to the data was generated in the same way, all the matrices have the same order

var.proteins<-proteins1

for(i in 2:ncol(var.proteins))
{
	var.proteins[,i]<-apply(cbind(proteins1[,i],proteins2[,i],proteins3[,i]),1,var)	
}
write.csv(var.proteins,file="variance noise data 0.5%.csv")

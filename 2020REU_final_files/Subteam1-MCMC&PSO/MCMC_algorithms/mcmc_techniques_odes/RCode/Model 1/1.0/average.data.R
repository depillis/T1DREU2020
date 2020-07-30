#source('average.data.R')

proteins1<-data.matrix(read.csv("noise data 1%_r1.csv", row.names=1))
proteins2<-data.matrix(read.csv("noise data 1%_r2.csv", row.names=1))
proteins3<-data.matrix(read.csv("noise data 1%_r3.csv", row.names=1))
## Due to the data was generated in the same way, all the matrices have the same order

proteins<-proteins1

proteins[,2:ncol(proteins)]<-(proteins1[,2:ncol(proteins1)]+proteins2[,2:ncol(proteins2)]+proteins3[,2:ncol(proteins3)])/3

write.csv(proteins,file="noise data 1%.csv")

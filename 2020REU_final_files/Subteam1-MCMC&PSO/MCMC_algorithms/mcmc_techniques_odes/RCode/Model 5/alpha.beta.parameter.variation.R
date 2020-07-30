observable.data<-data.matrix(read.csv(paste('Observable.data.csv',sep=''),sep=';'))

sd.matrix<-observable.data[,seq(4,ncol(observable.data),3)]*observable.data[,seq(5,ncol(observable.data),3)] #sd x factor, factor=1 data exits, factor=0 data does not exist
sd.vector<-c(sd.matrix)
sd.vector<-sd.vector[sd.vector!=0]
library(MASS)
sd.hyperparameter<-fitdistr(x=1/(sd.vector^2), densfun='gamma')
alpha.var<-sd.hyperparameter[[1]][1] #shape
 beta.var<-sd.hyperparameter[[1]][2]   #rate

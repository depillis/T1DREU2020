
############### Code
#With the collaboration of the team of Platform of Chemical and Biological Analysis using Computer Algebra Methods, who downloaded parameters values for more than 400 models.


directory.path<-'/.../Results_Biomodels/'
all.parameter<-c()
for(i in c(paste('00',1:9,sep=''),paste('0',10:99,sep=''),100:462)){
	temp.k <- tryCatch( read.csv(paste(directory.path,'/BIOMD0000000',i,'/Parameters_.txt',sep=''),sep='=',header=FALSE), error=function(e) 'error')
	# send error when there are reading problems, i.e. folder or file does not exit
	if(any(temp.k!='error'))all.parameter<-rbind(all.parameter, temp.k)
}
library(MASS)
k.ln <- fitdistr(x=log(na.exclude(all.parameter[all.parameter[,2]>0,2])) , densfun='normal')
k.ln.mean <- k.ln[[1]][1]
k.ln.sd <- k.ln[[1]][2]

prior.data<-c(k.ln.mean,k.ln.sd)
colnames(prior.data)<-c('MEAN','SD')
rownames(prior.data)<-c('k')
write.csv(prior.data,file='Hyperparameters values for kinetics parameter from Biomodels.csv')

pdf('histogram from biomodels.pdf')
x<-unique(sort( log(na.exclude(all.parameter[all.parameter[,2]>0,2])) ))
y<-(exp(-(unique(x)-k.ln.mean)^2/(2*pi^2)))/sqrt(2*pi*k.ln.sd^2) #density function
hist( log(na.exclude(all.parameter[all.parameter[,2]>0,2])), breaks=42 ,main=paste('Ln k ~ (',round(k.ln.mean,2),',',round(k.ln.sd,2),')',sep=''),ylab='Frequency',xlab='Ln k values',col='lightblue')
hist( log(na.exclude(all.parameter[all.parameter[,2]>0,2])), breaks=42 ,main=paste('Ln k ~ (',round(k.ln.mean,2),',',round(k.ln.sd,2),')',sep=''),ylab='Frequency',xlab='Ln k values',col='pink',freq=FALSE,ylim=c(0,1.1*max(y)))
lines(x,y)
dev.off()


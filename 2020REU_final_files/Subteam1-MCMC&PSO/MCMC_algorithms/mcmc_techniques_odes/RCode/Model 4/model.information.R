#source('model.information.R')
variance.vector<-(rep(5.6547672,18))^2  	#mass action
average.vector<-rep(-1.3697661,18)	#mass action
	
variance.vector[7]<-3.42346562^2  #k1f is michaelis menten kinetics, it is at 7th position
average.vector[7]<--2.59219616
###  New hyperparameter computation for Michaelis Mente kinetics, mean=-2.62 & sd=3.46


magnitude.variance<-sqrt(sum(variance.vector*variance.vector))
#print(str(average.vector))



##########reading the kinetic parameters##########
Ks<-data.matrix(read.csv(paste(path.array[k,1],n.l,'/set parameter solutions_egf path steady state daniel.csv',sep=''),sep=',',row.names=1))
Ks<-Ks[c('k1a','k1aBasic','k1b','k1c','k1d','k1e','k1f','k1g','k1r','k21','k22','k3','km2','km3','k_IRP_1Step','k_IRSiP_1Step','k_IRSiP_2Step','k_IRSiP_DosR'),]
random.parameter.vector<-names(Ks)
#Ks s a vector
fixed.par<-Ks
k.names<-c('k1a','k1aBasic','k1b','k1c','k1d','k1e','k1f','k1g','k1r','k21','k22','k3','km2','km3')
Ks.matrix<-t(array(Ks,dim=c(length(Ks),n.chains)))  #for metropolis
colnames(Ks.matrix) <- names(Ks)   ### for metropolis
log.current.sigma<-diag(0.1,length(Ks),length(Ks))  ### for adaptive
logik<-log(Ks)
###########################################################

################# reading experimental data #################
initial.values<-data.matrix(read.csv(paste(path.array[k,1],'/initial.species.values.csv',sep=''),sep=',',row.names=1))
observable.data<-data.matrix(read.csv(paste(path.array[k,1],'/Observable.data.csv',sep=''),sep=';'))
selected.species<-c('IR','IRins','IRp','IRiP','IRi','IRS','IRSiP','XX','XXp')
required.list<-list(initial.values,observable.data,selected.species)
#############################################################


######################################
## introduciendo los forcing values ##
######################################
times<-sort(unique(observable.data[,'time']))
forcing.matrix<-matrix(0,nrow=length(times),ncol=2)
forcing.matrix[,1]<-times
forcing.list<-list()
temp.list<-list()
insulin1.vector<-c(100,1.2,0,0.01,0.1,0.3,1,10,100)
for(z in 1:9){
	temp.list[[1]]<-forcing.matrix
	temp.list[[2]]<-forcing.matrix
	temp.list[[1]][,2]<-insulin1.vector[z]
	if(z==2) temp.list[[2]][temp.list[[2]][,1]>=4,2]<-8.8
	forcing.list[[z]]<-temp.list
}

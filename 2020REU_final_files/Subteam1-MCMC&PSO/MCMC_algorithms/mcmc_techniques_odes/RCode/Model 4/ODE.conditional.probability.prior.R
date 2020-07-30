#ode.conditional.probability<-function(ks,required.data,alpha.var,beta.var,selected.species,forcing,dllname)
ode.conditional.probability<-function(ks,k.names,hyperparameter,required.data,forcing.list,dllname)
# hyperparameter --> alpha.var, beta.var
# required.data --> initial.data, observable.data, selected.species
#ode.conditional.probability<-function(ks,k.names,alpha.var,beta.var,initial.data,observable.data,selected.species,forcing.list,dllname)
{
library(deSolve)
alpha.var<-hyperparameter[[1]]
 beta.var<-hyperparameter[[2]]
initial.data<-required.data[[1]]
observable.data<-required.data[[2]]
selected.species<-required.data[[3]]
parms<-ks[k.names]
#print(parms)
xstart<-initial.data[selected.species,]  #cambiado del 22 de julio de  2019
#print(xstart)
times<-sort(unique(observable.data[,'time']))

k_IRP_1Step<-ks['k_IRP_1Step']   		#scaling parameter
k_IRSiP_1Step<-ks['k_IRSiP_1Step']		#scaling parameter
k_IRSiP_2Step<-ks['k_IRSiP_2Step']		#scaling parameter
k_IRSiP_DosR <-ks['k_IRSiP_DosR']		#scaling parameter
#print(dllname)   # <<<<<<<<<----------- for testing

simulated.data<-c()
exponential<-0
#for(z in 1:2){   # <<<<<<<<<----------- for testing
for(z in 1:9){
	forcing<-forcing.list[[z]]
	time.obs<-data.matrix( observable.data[observable.data[,1]==z,'time'])
	#print(dllname)   # <<<<<<<<<----------- for testing
	#print(str(ode(xstart[selected.species],times,func = "derivs",parms = parms, dllname = dllname, initforc="forcc",forcings = forcing, initfunc = "initmod", nout = 1,outnames = "Sum")))
	out1<-tryCatch(ode(xstart[selected.species],times,func = "derivs",parms = parms, dllname = dllname, initforc="forcc",forcings = forcing, initfunc = "initmod", nout = 1,outnames = "Sum"),error=function(e) -Inf)
	#print(str(out1))   # <<<<<<<<<----------- for testing
	#print(paste('----------------',z,'-------------------'))   # <<<<<<<<<----------- for testing
	if(length(out1)!=1){
		if(!all(times %in% out1[,1])){	exponential<--1e+308
		}else{
			temp.obs<-as.matrix(out1[out1[,'time']%in%time.obs,c('time','IRp','IRiP','IRSiP')])
			if(ncol(temp.obs)==1)temp.obs<-t(temp.obs)

			#IRS1_P<-k_IRSiP_DosR*temp.obs[,'IRSiP']/sigmaYDosR
			IRS1_P<-k_IRSiP_DosR*temp.obs[,'IRSiP']
			IR1_P<-rep(0,nrow(temp.obs))
			if(z==1){
				#IR1_P <-k_IRP_1Step*(temp.obs[,'IRiP'] + temp.obs[,'IRp'])/sigmaY1TimR
				IR1_P <-k_IRP_1Step*(temp.obs[,'IRiP'] + temp.obs[,'IRp'])
				#IRS1_P<-k_IRSiP_1Step*temp.obs[,'IRSiP']/sigmaY2TimR
				IRS1_P<-k_IRSiP_1Step*temp.obs[,'IRSiP']
			}
			if(z==2)IRS1_P<-k_IRSiP_2Step*temp.obs[,'IRSiP']

			temp.obs2<-cbind(z,temp.obs[,'time'],IR1_P,IRS1_P)
		
			if(ncol(temp.obs2)==1){names(temp.obs2)[1:2]<-c('Model','time')
			}else{colnames(temp.obs2)[1:2]<-c('Model','time')}

			simulated.data<-rbind(simulated.data,temp.obs2)
		}
	}else{exponential<--1e+308}
	#print(tail(simulated.data[,1:3]))   # <<<<<<<<<----------- for testing

}

rownames(simulated.data)<-NULL
#print(tail( simulated.data[,c('IR1_P','IRS1_P')]))   # <<<<<<<<<----------- for testing
#print(tail(observable.data[,c('IR1_P','IRS1_P')]))   # <<<<<<<<<----------- for testing
if(exponential==0) {
	#IR1_P	IR1_P_std	IR1_P_ExpError	IRS1_P	IRS1_P_std	IRS1_P_ExpError

	#error.matrix<-sum(((observable.data[,c('IR1_P','IRS1_P')]-simulated.data[,c('IR1_P','IRS1_P')])^2)/(observable.data[,c('IR1_P_std','IRS1_P_std')])^2)
	#exponential<- -1*sum(error.matrix)

	sqrt.dif<-(observable.data[,c('IR1_P','IRS1_P')]-simulated.data[,c('IR1_P','IRS1_P')])^2
	#sqrt.dif<-(observable.data[1:24,c('IR1_P','IRS1_P')]-simulated.data[,c('IR1_P','IRS1_P')])^2 ## porque solo esta tratamiento 1 y 2
	error.matrix<-(2*alpha.var+1/2)*log(1+sqrt.dif/(2*beta.var))			#22.09.2014
	#print(str(error.matrix))    # <<<<<<<<<----------- for testing
	#print(quantile(error.matrix))   # <<<<<<<<<----------- for testing

	error.matrix[error.matrix=="NaN"]<-Inf
	exponential<--sum(error.matrix)

}
########################
##### Observables ######
#### Model 1 ####
# IR1_P<-k_IRP_1Step*(IRiP + IRp)
# IRS1_P<-IRSiP*k_IRSiP_1Step
#### Model 2 ####
# IRS1_P<-IRSiP*k_IRSiP_2Step
#### Model 3 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 4 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 5 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 6 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 7 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 8 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 9 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#################
# parámetros de escalamiento
#define k_IRP_1Step parms[14]
#define k_IRSiP_1Step parms[15]
#define k_IRSiP_2Step parms[16]
#define k_IRSiP_DosR parms[17]
#define sigmaY1TimR parms[18]
#define sigmaY2Step parms[19]
#define sigmaY2TimR parms[20]
#define sigmaYDosR parms[21]

return(exponential)
      
}

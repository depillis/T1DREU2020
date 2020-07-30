#source("MCMC_metropolis_likelihood trial_experimental data.R")
n.chains=1
likelihood.choice.vector<-c('','.varpriors','.kpriors','.k_varpriors')
#path.root='E:/'
#path.root='/R codes model/'
#path.root<-'C:/R codes model/'

library('snow')
library('parallel')
library(coda)
library('adaptMCMC')
library('MASS')
#library(mcmc)

likelihood.choice<-likelihood.choice.vector[2]  			## selecting priors for the objetctive function
n.set.vector<-c(25,25,13,18,23)								## amount of parameter to be computed per model
n.part<-paste('single.adaptive',likelihood.choice,sep='')

k<-1    #k=1:3, for the five models								## variable for selected model
n.set<-n.set.vector[k]
load('variance.hyperparameter.RData')		 ## parameter.list		#22.09.2014

path.array<-array(0,dim=c(5,5))

path.array[1,1]<-paste(path.root,'Model 1/',sep='')
path.array[2,1]<-paste(path.root,'Model 2/',sep='')
path.array[3,1]<-paste(path.root,'Model 3/',sep='')
path.array[4,1]<-paste(path.root,'Model 4/',sep='')
path.array[5,1]<-paste(path.root,'Model 5',sep='')

path.array[,3]<-paste('ODE.conditional.probability.R',sep='')
if(likelihood.choice=='.varpriors'|likelihood.choice=='.k_varpriors')path.array[,3]<-paste('ODE.conditional.probability.prior.R',sep='')

path.array[1,4]<-'mapkk_feedback'
path.array[2,4]<-'jak'
path.array[3,4]<-'smad'
path.array[4,4]<-'Brannmark'
path.array[5,4]<-'Raia2011'

setwd(path.array[k,1])
print(c(path.array[k,1],path.array[k,2]))
system(paste('R CMD SHLIB ',path.array[k,4],'.c',sep='')) #ozi
dyn.load(paste(path.array[k,4],".so",sep='')) #ozi

################### variation information ####################
load('variance.hyperparameter.RData')		 ## parameter.list		#22.09.2014
##########################################################


### calling.function_a is using for computing the log likelihood
if(likelihood.choice==''){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		#names(par)<-logik.names
		current.proba.exp<-ode.conditional.probability(ks=par,required.data=all.data)
		print(c(paste('i: ',i,sep=''),current.proba.exp))
		return(current.proba.exp)
	}
}
if(likelihood.choice=='.kpriors'){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		#names(par)<-logik.names
		current.proba.exp<-ode.conditional.probability(ks=par,required.data=all.data)
		proba.exp.prior<-current.proba.exp-(n.set/2)*log(2*pi)-1/2*log(magnitude.variance)-1/2*sum(((log(par)-average.vector)^2)/variance.vector) #i am not taking into account n, cause there are not pre information
		print(c(paste('i: ',i,sep=''),current.proba.exp,proba.exp.prior))
		return(proba.exp.prior)
	}
}
if(likelihood.choice=='.varpriors'){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		#names(par)<-logik.names
		current.proba.exp<-ode.conditional.probability(ks=par,required.data=all.data,alpha.var=alpha.var,beta.var=beta.var)
		print(c(paste('i: ',i,sep=''),current.proba.exp))
		return(current.proba.exp)
	}
}
if(likelihood.choice=='.k_varpriors'){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		#names(par)<-logik.names
		current.proba.exp<-ode.conditional.probability(ks=par,required.data=all.data,alpha.var=alpha.var,beta.var=beta.var)
		proba.exp.prior<-current.proba.exp-(n.set/2)*log(2*pi)-1/2*log(magnitude.variance)-1/2*sum(((log(par)-average.vector)^2)/variance.vector) #i am not taking into account n, cause there are not pre information
		print(c(paste('i: ',i,sep=''),current.proba.exp,proba.exp.prior))
		return(proba.exp.prior)
	}
}
##################################################
##################################################
variance.list<-list()
average.list<-list()

variance.list[[1]]<-c(3.42346562,log(25),4.41506944,rep(3.42346562,19)))^2
average.list[[1]]<-c(-2.59219616,log(2),-7.21098853,rep(-2.59219616,19))
variance.list[[2]]<-(rep(5.6547672,25))^2
average.list[[2]]<-rep(-1.3697661,25)
variance.list[[3]]<-(rep(3.42346562,13))^2
average.list[[3]]<-rep(-2.59219616,13)

variance.vector<-variance.list[[k]]
average.vector<-average.list[[k]]
magnitude.variance<-sqrt(sum(variance.vector*variance.vector))

for(n.l in c('5.0','10.0','15.0','25.0')){
	i<-1
	path.array[,2]<-n.l
	n.l2<-''
	if(n.l=='5.0'){n.l2<-'d'}
	if(n.l=='10.0'){n.l2<-'e'}
	if(n.l=='15.0'){n.l2<-'f'}
	if(n.l=='25.0'){n.l2<-'g'}
	path.array[,5]<-n.l2

	source(paste(path.array[k,1],path.array[k,3],sep=''))  # ODE.probability...R

	################  Reading model information ###########
	load(paste(path.array[k,1],'data/input.data.',path.array[k,5],'.RData',sep=''))
	k.data<-all.data[[6]]  #kinetic parameter computed in steady state
	if(n.l=='10.0'&k==2)k.data<-all.data[[7]]
	k.id<-colnames(k.data)
	Ks<-k.data[,1]
	log.current.sigma<-diag(0.1,length(Ks),length(Ks))  ### for adaptive
	logik<-log(Ks)
	###########################################################

	################### variation information ####################
	alpha.var<-as.numeric(parameter.list[[n.l2]]['fit.parameters',1])    #alpha col1			#08.10.2014
	beta.var<-as.numeric(parameter.list[[n.l2]]['fit.parameters',2]) 	 #alpha col2			#08.10.2014
	##########################################################

	ptm<-proc.time()
	all.samp<-MCMC(calling.function_a, n=1001000, init=logik,n.start=1000,scale=log.current.sigma,adapt=TRUE, acc.rate=0.234)
	ptm2<-proc.time()-ptm

	save(all.samp,ptm2,file=paste(path.array[k,1],path.array[k,2],'/all_',n.part,'.',n.chains,' chains_output.RData',sep=''))
	samp.added.data<-list(exp(all.samp$samples),all.samp$log.p)
	save(samp.added.data,ptm2,file=paste(path.array[k,1],path.array[k,2],'/all_',n.part,'.',n.chains,' chains_samples.RData',sep=''))

}


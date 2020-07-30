#source("MCMC_adaptive_likelihood trial_experimental data.R")
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

likelihood.choice<-likelihood.choice.vector[1]  			## selecting priors for the objetctive function
n.set.vector<-c(25,25,13,18,23)								## amount of parameter to be computed per model
n.part<-paste('single.adaptive',likelihood.choice,sep='')	## text for the RData file name

n.l<-'realistic'											## variable for noise level folder, here it is just experimental daata
k<-4														## variable for selected model
n.set<-n.set.vector[k]

path.array<-array(0,dim=c(5,5))

path.array[1,1]<-paste(path.root,'Model 1/',sep='')
path.array[2,1]<-paste(path.root,'Model 2/',sep='')
path.array[3,1]<-paste(path.root,'Model 3/',sep='')
path.array[4,1]<-paste(path.root,'Model 4/',sep='')
path.array[5,1]<-paste(path.root,'Model 5/',sep='')

path.array[,3]<-paste('ODE.conditional.probability.R',sep='')
if(likelihood.choice=='.varpriors'|likelihood.choice=='.k_varpriors')path.array[,3]<-paste('ODE.conditional.probability.prior.R',sep='')

path.array[1,4]<-'mapkk_feedback'
path.array[2,4]<-'jak'
path.array[3,4]<-'smad'
path.array[4,4]<-'Brannmark'
path.array[5,4]<-'Raia2011'

n.l2<-''

if(n.l=='5.0'){n.l2<-'d'}
if(n.l=='10.0'){n.l2<-'e'}
if(n.l=='15.0'){n.l2<-'f'}
if(n.l=='25.0'){n.l2<-'g'}

path.array[,5]<-n.l2

setwd(path.array[k,1])
source('model.information.R')  ### read all the variables to be used


################### variation information ####################
## Hyperparameters are computed from the variance data given in the repository made by Hass et. all
## Hass H, Loos C, Alvarez ER, Timmer J, Hasenauer J, Kreutz C. Benchmark problems for dynamic modeling of intracellular processes. BioRxiv. 2018;: p. 404590.
source('variance computation.R')
hyperparameter<-list(alpha.var,beta.var)
##########################################################


### calling.function_a is using for computing the log likelihood ###
if(likelihood.choice==''){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		current.proba.exp<-ode.conditional.probability(ks=par,k.names=k.names,initial.data=initial.values,observable.data,selected.species,forcing.list,dllname=path.array[k,4])
		print(c(paste('i: ',i,sep=''),current.proba.exp))
		return(current.proba.exp)
	}
}
if(likelihood.choice=='.kpriors'){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		current.proba.exp<-ode.conditional.probability(ks=par,k.names=k.names,initial.data=initial.values,observable.data,selected.species,forcing.list,dllname=path.array[k,4])
		proba.exp.prior<-current.proba.exp-(n.set/2)*log(2*pi)-1/2*log(magnitude.variance)-1/2*sum(((log(par)-average.vector)^2)/variance.vector) #i am not taking into account n, cause there are not pre information
		print(c(paste('i: ',i,sep=''),current.proba.exp,proba.exp.prior))
		return(proba.exp.prior)
	}
}
if(likelihood.choice=='.varpriors'){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		current.proba.exp<-ode.conditional.probability(ks=par,k.names=k.names,hyperparameter=hyperparameter,required.data=required.list,forcing.list,dllname=path.array[k,4])
		print(c(paste('i: ',i,sep=''),current.proba.exp))
		return(current.proba.exp)
	}
}
if(likelihood.choice=='.k_varpriors'){
	calling.function_a<-function(state)
	{
		par=exp(state) #for tempering is exp(state[-1])
		current.proba.exp<-ode.conditional.probability(ks=par,k.names=k.names,hyperparameter=hyperparameter,required.data=required.list,forcing.list,dllname=path.array[k,4])
		proba.exp.prior<-current.proba.exp-(n.set/2)*log(2*pi)-1/2*log(magnitude.variance)-1/2*sum(((log(par)-average.vector)^2)/variance.vector) #i am not taking into account n, cause there are not pre information
		print(c(paste('i: ',i,sep=''),current.proba.exp,proba.exp.prior))
		return(proba.exp.prior)
	}
}
######################


i<-1
print(c(path.array[k,1],path.array[k,2]))
system(paste('R CMD SHLIB ',path.array[k,4],'.c',sep='')) #ozi
dyn.load(paste(path.array[k,4],".so",sep='')) #ozi
source(paste(path.array[k,1],path.array[k,3],sep=''))  # ODE.probability...R


ptm<-proc.time()
all.samp<-MCMC(calling.function_a, n=1001000, init=logik,n.start=1000,scale=log.current.sigma,adapt=TRUE, acc.rate=0.234)
ptm2<-proc.time()-ptm

print(ptm2)

save(all.samp,ptm2,file=paste( path.array[k,1],'/realistic/all_',n.part,'.',n.chains,' chains_output.RData',sep=''))
samp.added.data<-list(exp(all.samp$samples),all.samp$log.p)
save(samp.added.data,ptm2,file=paste( path.array[k,1],'/realistic/all_',n.part,'.',n.chains,' chains_samples.RData',sep=''))





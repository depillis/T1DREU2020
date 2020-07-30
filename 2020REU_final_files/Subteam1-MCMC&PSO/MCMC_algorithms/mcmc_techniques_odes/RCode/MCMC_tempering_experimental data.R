#source('MCMC_tempering_experimental data.R')
#path.root='E:/'
#path.root='/R codes model/'
#path.root<-'C:/R codes model/'
n.part<-'tempering'

#library('adaptMCMC')
library(mcmc)
library('snow')
library('parallel')
library(coda)
library(MASS)
n.set.vector<-c(25,25,13,18,23)	

##### calling.function_a is using for computing the log likelihood
calling.function<-function(state)
{
	par=exp(state[-1]) #for tempering is exp(state[-1])
	names(par)<-logik.names
	current.proba.exp<-ode.conditional.probability(ks=par,k.names=k.names,initial.data=initial.values,observable.data,selected.species,forcing.list,dllname=path.array[k,4])
	print(current.proba.exp)
	return(current.proba.exp)
}
##########################

path.array<-array(0,dim=c(5,5))

path.array[1,1]<-paste(path.root,'Model 1/',sep='')
path.array[2,1]<-paste(path.root,'Model 2/',sep='')
path.array[3,1]<-paste(path.root,'Model 3/',sep='')
path.array[4,1]<-paste(path.root,'Model 4/',sep='')
path.array[5,1]<-paste(path.root,'Model 5/',sep='')

path.array[,3]<-paste('ODE.conditional.probability.R',sep='')

path.array[1,4]<-'mapkk_feedback'
path.array[2,4]<-'jak'
path.array[3,4]<-'smad'
path.array[4,4]<-'Brannmark'
path.array[5,4]<-'Raia2011'

for(n.chains in c(5,10,15)){#'5.0',10.0,'15.0','25.0'
	n.l<-'realistic'  #there are not noise levels

	path.array[,2]<-n.l
	path.array[,3]<-'ODE.conditional.probability.R'

	n.l2<-''

	if(n.l=='0.1'){n.l2<-'a'}
	if(n.l=='0.5'){n.l2<-'b'}
	if(n.l=='1.0'){n.l2<-'c'}
	if(n.l=='5.0'){n.l2<-'d'}
	if(n.l=='10.0'){n.l2<-'e'}
	if(n.l=='15.0'){n.l2<-'f'}
	if(n.l=='25.0'){n.l2<-'g'}
	path.array[,5]<-n.l2
	k<-4
	setwd(path.array[k,1])
	source('model.information.R')  ### read all the variables to be used

	################### variation information ####################
	## Hyperparameters are computed from the variance data given in the repository made by Hass et. all
	## Hass H, Loos C, Alvarez ER, Timmer J, Hasenauer J, Kreutz C. Benchmark problems for dynamic modeling of intracellular processes. BioRxiv. 2018;: p. 404590.
	source('variance computation.R')
	hyperparameter<-list(alpha.var,beta.var)
	##########################################################

	print(c(path.array[k,1],path.array[k,2]))
	system(paste('R CMD SHLIB ',path.array[k,4],'.c',sep='')) #ozi
	dyn.load(paste(path.array[k,4],".so",sep='')) #ozi
	source(paste(path.array[k,1],path.array[k,3],sep=''))  # ODE.probability...R
	
	#############generate the neighbors matrix###################
	neighbors <- matrix(TRUE,n.chains,n.chains)
	###########################################################

	logik.names<-names(logik)
	ptm<-proc.time()
	all.samp<-temper(calling.function, initial = log(Ks.matrix), neighbors = neighbors, nbatch = 10010,blen = 100, nspac = 1, scale = 0.01, parallel = TRUE, debug = TRUE)
	ptm2<-proc.time()-ptm
	print(ptm2)

	save(all.samp,ptm2,file=paste(path.array[k,1],path.array[k,2],'/all_',n.part,'.',n.chains,' chains_output.RData',sep=''))
	chains.together.data<-exp(all.samp$state)
	save(chains.together.data,ptm2,file=paste(path.array[k,1],path.array[k,2],'/all_',n.part,'.',n.chains,' chains_samples.RData',sep=''))  ## this part i should do later
}


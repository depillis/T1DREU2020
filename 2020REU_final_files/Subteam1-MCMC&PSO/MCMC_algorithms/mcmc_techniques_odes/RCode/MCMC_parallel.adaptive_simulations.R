#source('MCMC_parallel.adaptive_experimental data.R')
#path.root='E:/'
#path.root='/R codes model/'
#path.root<-'C:/R codes model/'
n.part<-'parallel.adaptive'

library('snow')
library('parallel')
library(coda)
library('adaptMCMC')
library('MASS')
n.set.vector<-c(25,25,13,18,23)	

##### calling.function_a is using for computing the log likelihood
calling.function<-function(state)
{
	par=exp(state) #for tempering is exp(state[-1])
	current.proba.exp<-ode.conditional.probability(ks=par,required.data=all.data)
	print(current.proba.exp)
	return(current.proba.exp)
}
##################################################
#################################################

######  Adaptive parallelization function###
## It is used when there are problem with allocate the memory
samp.fun0<-function(x)
{
	samp.initial<-MCMC(calling.function, n=1001000, init=logik,n.start=1000,scale=log.current.sigma,adapt=TRUE, acc.rate=0.234)
	return(samp.initial)
}
####################################


path.array<-array(0,dim=c(5,5))

path.array[1,1]<-paste(path.root,'20130215_noise analysis/',sep='')
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

k<-1  #models 1:3
setwd(path.array[k,1])

print(c(path.array[k,1],path.array[k,2]))
system(paste('R CMD SHLIB ',path.array[k,4],'.c',sep='')) #ozi
dyn.load(paste(path.array[k,4],".so",sep='')) #ozi

for(n.chains in c(1,5,10,15)){
	#n.chains<-15
	for(n.l in c('0.1','0.5','1.0','5.0','10.0','15.0','25.0')){

		path.array[,2]<-n.l
		n.l2<-''

		if(n.l=='0.1'){n.l2<-'a'}
		if(n.l=='0.5'){n.l2<-'b'}
		if(n.l=='1.0'){n.l2<-'c'}
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


		ptm<-proc.time()
		samp.1<-lapply(1:n.chains, function(x) mcparallel(samp.fun0(x),name=x))
		all.samp<-mccollect(samp.1)
		ptm2<-proc.time() - ptm

		##If there are not memory allocation problems it can be used 
		## MCMC.parallel(calling.function, n=1001000, init=logik,n.start=1000,n.chain=n.chains,n.cpu=n.chains, scale=log.current.sigma,adapt=TRUE, acc.rate=0.234)  ### just watn to see is for time n=251000
		samp.added.data<-lapply(1:n.chains, function(x) list(exp(all.samp[[x]]$samples),all.samp[[x]]$log.p))
		if(n.chains==1){all.samp<-all.samp[[1]]}
		save(all.samp,ptm2,file=paste(path.array[k,1],path.array[k,2],'/all_',n.part,'.',n.chains,' chains_output.RData',sep=''))
		save(samp.added.data,ptm2,file=paste(path.array[k,1],path.array[k,2],'/all_',n.part,'.',n.chains,' chains_samples.RData',sep=''))
	}
}


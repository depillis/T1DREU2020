#setwd("/R codes model/Model 3")
#source("alpha.beta.parameter.variation.R")
#k<-3

library(mcmc)
parameter.list<-list()
opt.parameter.list<-list()
mcmc.parameter.list<-list()
fit.parameter.list<-list()
parameter.list[['a']]<-NULL
parameter.list[['b']]<-NULL
parameter.list[['c']]<-NULL
parameter.list[['d']]<-NULL
parameter.list[['e']]<-NULL
parameter.list[['f']]<-NULL
parameter.list[['g']]<-NULL

selected.species.list[[1]]<-c("MKKK","MKKK_P","MKK","MKK_P","MKK_PP","MAPK","MAPK_P","MAPK_PP")
selected.species.list[[2]]<-c('JAK2','pJAK2','EpoR','pEpoR','SHP1','mSHP1','actSHP1','SOS','mSOS','Raf','pRaf','MEK2','pMEK2','MEK1','pMEK1','ppMEK2','ppMEK1','ERK1','pERK1','ERK2','pERK2','ppERK1','ppERK2','pSOS')
selected.species.list[[3]]<-c('receptors','R_Smad_cyt','R_Smad_P_cyt','Smad4_cyt','R_Smad_P_Smad4_cyt','R_Smad_P_Smad4_nuc','R_Smad_nuc','R_Smad_P_nuc','Smad4_nuc','Pi')


pdf('Inverse Sigma Data Distribution - 15.0 noise.pdf',width=18,height=7)

for(adapt.choice in c('a','b','c','d','e','f','g')){
	load(paste(getwd(),'/data/input.data.',adapt.choice,'.RData',sep=''))

	initial.data<-all.data[[1]]
	data.matrix1<-all.data[[2]]
	data.matrix2<-all.data[[3]]
	data.matrix3<-all.data[[4]]
	sigma.data<-all.data[[5]]   ### variance of noise level

	selected.species<-selected.species.list[[k]]
	inverse.sigma.data<-1/(sigma.data[,selected.species])

	min.function<-function(x)
	{
		alpha=x[1]    	##  JAK2_phosphorylation_by_Epo   		##
		beta=x[2]    	##  EpoR_phosphorylation_by_pJAK2 		##
		sum((inverse.sigma.data-alpha/beta)^2)+sum(((inverse.sigma.data-alpha/beta)^2-(alpha/beta^2))^2)
	}

	min.function2<-function(x)
	{
		alpha=exp(x[1])    	##  JAK2_phosphorylation_by_Epo   		##
		beta=exp(x[2])    	##  EpoR_phosphorylation_by_pJAK2 		##
		return(-(sum((inverse.sigma.data-alpha/beta)^2)+sum(((inverse.sigma.data-alpha/beta)^2-(alpha/beta^2))^2)))
	}


	ui.m<-diag(1,2)
	ci.v<-rep(1e-7,2)
	P<-constrOptim(rep(10,2),min.function,NULL,ui=ui.m,ci=ci.v)

	parameters<-P$par

	out <- metrop(min.function2, log(c(10,1)), 2000,debug=TRUE)
	mcmc.parameters<-exp(apply(out$batch[1001:2000,],2,mean))
	
	out <- metrop(min.function2, log(c(10,1)), 500000,debug=TRUE)
	mcmc2.parameters<-exp(apply(out$batch[250001:500000,],2,mean))

	real.distribution<-c(inverse.sigma.data)
	distribution1<-rgamma(length(real.distribution), shape=parameters[1], rate = parameters[2])
	distribution2<-rgamma(length(real.distribution), shape=mcmc.parameters[1], rate = mcmc.parameters[2])
	distribution2.2<-rgamma(length(real.distribution), shape=mcmc2.parameters[1], rate = mcmc2.parameters[2])

	library(MASS)
	fit.parameters<-tryCatch(fitdistr(real.distribution,'gamma'), error = function(e) c('error','error'))
	distribution3<-tryCatch(rgamma(length(real.distribution), shape=fit.parameters[[1]][1], rate = fit.parameters[[1]][2]),error = function(e) 0)

	par(mfrow=c(1,5)) #las=2 perpendicular
	hist(distribution1,main=paste('Optimization',adapt.choice,sep=' - '))
	hist(distribution2,main=paste('Metropolis Hastings',adapt.choice,sep=' - '))
	hist(distribution2.2,main=paste('Metropolis Hastings +iter',adapt.choice,sep=' - '))
	hist(distribution3,main=paste('Fit distribution',adapt.choice,sep=' - '))
	hist(real.distribution,main=paste('Inverse Sigma Data',adapt.choice,sep=' - '))

	hist(distribution1,main=paste('Optimization',adapt.choice,sep=' - '),xlim=range(real.distribution))
	hist(distribution2,main=paste('Metropolis Hastings',adapt.choice,sep=' - '),xlim=range(real.distribution))
	hist(distribution2.2,main=paste('Metropolis Hastings +iter',adapt.choice,sep=' - '),xlim=range(real.distribution))
	hist(distribution3,main=paste('Fit Distribution',adapt.choice,sep=' - '),xlim=range(real.distribution))
	hist(real.distribution,main=paste('Inverse Sigma Data',adapt.choice,sep=' - '),xlim=range(real.distribution))

	alpha=(parameters[1]+mcmc.parameters[1])/2
	beta=(parameters[2]+mcmc.parameters[2])/2

	parameter.list[[adapt.choice]]<-c(alpha,beta)
	fit.parameters<-c(alpha=fit.parameters[[1]][1],beta=fit.parameters[[1]][2])
	mcmc.parameters<-c(alpha=mcmc.parameters[1],beta=mcmc.parameters[2])
	opt.parameters<-c(alpha=parameters[1],beta=parameters[2])
	mcmc2.parameters<-c(alpha=mcmc2.parameters[1],beta=mcmc2.parameters[2])
	parameter.list[[adapt.choice]]<-rbind(fit.parameters,opt.parameters, mcmc.parameters,mcmc2.parameters)
}

dev.off()
save(parameter.list, file='variance.hyperparameter.RData')

#setwd("/R codes model/Model 2")
#source("alpha.beta.parameter.variation.R")

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
pdf('Inverse Sigma Data Distribution - 15.0 noise.pdf',width=18,height=7)
#for(adapt.choice in c('d','e','f','g')){
for(adapt.choice in c('a','b','c','d','e','f','g')){

load(paste('/home/abidata/Gloria/2011-2012/20130225_noise analysis2/data/input.data.',adapt.choice,'.RData',sep=''))
#load(paste('/home/abidata/Gloria/2011-2012/20130215_noise analysis/data/input.data.',adapt.choice,'.RData',sep=''))

initial.data<-all.data[[1]]
data.matrix1<-all.data[[2]]
data.matrix2<-all.data[[3]]
data.matrix3<-all.data[[4]]
sigma.data<-all.data[[5]]

selected.species<-c('JAK2','pJAK2','EpoR','pEpoR','SHP1','mSHP1','actSHP1','SOS','mSOS','Raf','pRaf','MEK2','pMEK2','MEK1','pMEK1','ppMEK2','ppMEK1','ERK1','pERK1','ERK2','pERK2','ppERK1','ppERK2','pSOS')

inverse.sigma.data<-1/(sigma.data[,selected.species])


min.function<-function(x)
{
#print(x)
alpha=x[1]    	##  JAK2_phosphorylation_by_Epo   		##
beta=x[2]    	##  EpoR_phosphorylation_by_pJAK2 		##

sum((inverse.sigma.data-alpha/beta)^2)+sum(((inverse.sigma.data-alpha/beta)^2-(alpha/beta^2))^2)
}

min.function2<-function(x)
{
#print(exp(x))
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
#out <- metrop(min.function2, log(c(10,1)), 500000,debug=TRUE)
#mcmc2.parameters<-exp(apply(out$batch[250001:500000,],2,mean))
out <- metrop(min.function2, log(c(10,1)), 2000000,debug=TRUE)
mcmc2.parameters<-exp(apply(out$batch[1000001:2000000,],2,mean))

real.distribution<-c(inverse.sigma.data)
distribution1<-rgamma(length(real.distribution), shape=parameters[1], rate = parameters[2])
distribution2<-rgamma(length(real.distribution), shape=mcmc.parameters[1], rate = mcmc.parameters[2])
distribution2.2<-rgamma(length(real.distribution), shape=mcmc2.parameters[1], rate = mcmc2.parameters[2])


library(MASS)
#fit.parameters<-tryCatch(fitdistr(real.distribution,'gamma'), error = function(e) c('error','error'))
fit.parameters<-fitdistr(real.distribution,'gamma',lower=0.0011)
print(fit.parameters)
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


pdf('Inverse Sigma Data Distribution - 15.0 noise.pdf',width=18,height=7)
for(adapt.choice in c('a','b','c','d','e','f','g')){
load(paste('/home/abidata/Gloria/2011-2012/20130225_noise analysis2/data/input.data.',adapt.choice,'.RData',sep=''))
sigma.data<-all.data[[5]]
selected.species<-c('JAK2','pJAK2','EpoR','pEpoR','SHP1','mSHP1','actSHP1','SOS','mSOS','Raf','pRaf','MEK2','pMEK2','MEK1','pMEK1','ppMEK2','ppMEK1','ERK1','pERK1','ERK2','pERK2','ppERK1','ppERK2','pSOS')
real.distribution<-c(1/(sigma.data[,selected.species]))
distribution1<-rgamma(length(real.distribution), shape=as.numeric(parameter.list[[adapt.choice]][2,1]),rate=as.numeric(parameter.list[[adapt.choice]][2,2]))
distribution2<-rgamma(length(real.distribution), shape=as.numeric(parameter.list[[adapt.choice]][3,1]),rate=as.numeric(parameter.list[[adapt.choice]][3,2]))
distribution2.2<-rgamma(length(real.distribution), shape=as.numeric(parameter.list[[adapt.choice]][4,1]),rate=as.numeric(parameter.list[[adapt.choice]][4,2]))
distribution3<-rgamma(length(real.distribution), shape=as.numeric(parameter.list[[adapt.choice]][1,1]),rate=as.numeric(parameter.list[[adapt.choice]][1,2]))



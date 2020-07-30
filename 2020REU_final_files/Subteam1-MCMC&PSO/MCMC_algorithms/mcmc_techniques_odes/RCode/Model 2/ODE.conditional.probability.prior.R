ode.conditional.probability<-function(ks,required.data,alpha.var,beta.var)
{
library(deSolve)

initial.data<-required.data[[1]]
data.matrix1<-required.data[[2]]
data.matrix2<-required.data[[3]]
data.matrix3<-required.data[[4]]
sigma.data<-required.data[[5]]


xstart<-c(initial.data['JAK2'],initial.data['pJAK2'],initial.data['Epo'],initial.data['EpoR'],initial.data['pEpoR'],initial.data['SHP1'],initial.data['mSHP1'],initial.data['Delay01_mSHP1'],initial.data['Delay02_mSHP1'],initial.data['Delay03_mSHP1'],initial.data['Delay04_mSHP1'],initial.data['Delay05_mSHP1'],initial.data['Delay06_mSHP1'],initial.data['Delay07_mSHP1'],initial.data['Delay08_mSHP1'],initial.data['actSHP1'],initial.data['SOS'],initial.data['mSOS'],initial.data['Raf'],initial.data['pRaf'],initial.data['MEK2'],initial.data['pMEK2'],initial.data['MEK1'],initial.data['pMEK1'],initial.data['ppMEK2'],initial.data['ppMEK1'],initial.data['ERK1'],initial.data['pERK1'],initial.data['ERK2'],initial.data['pERK2'],initial.data['ppERK1'],initial.data['ppERK2'],initial.data['pSOS'])
##compartment=4e-12

parms<-c(k1=ks[1],k2=ks[2],k3=ks[3],k4=ks[4],k13=ks[5],k14=ks[6],k15=ks[7],k16=ks[8],k17=ks[9],k18=ks[10],k19=ks[11],k20=ks[12],k21=ks[13],k22=ks[14],k23=ks[15],k24=ks[16],k26=ks[17],k28=ks[18],k29=ks[19],k32=ks[20],k33=ks[21],k36=ks[22],k38=ks[23],k40=ks[24],k42=ks[25],cell=1)

#print(parms)

times<-c(0,0.25,0.5,1,1.5,3,4.5,6,7.5,9,10.5,12,18,24,36,48,100)

out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "jak", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes"),error=function(e) -Inf) 

out1<-data.frame(out1)
cols<-ncol(out1)

sim.length<-nrow(out1)
tim.length<-length(times)
if(sim.length!=tim.length){
	exponential<--1e+300  # change it to Inf on february 9, 2019 
	#exponential<--Inf
}
else if(!all(times %in% out1[,1])){
	exponential<--1e+300
	#exponential<--Inf
}else{
	#print('computing the likelihood')
	selected.species<-c('JAK2','pJAK2','EpoR','pEpoR','SHP1','mSHP1','actSHP1','SOS','mSOS','Raf','pRaf','MEK2','pMEK2','MEK1','pMEK1','ppMEK2','ppMEK1','ERK1','pERK1','ERK2','pERK2','ppERK1','ppERK2','pSOS')

	#------------------------------ Computation with data points ----------------------------------------
	sqrt.dif<-((out1[,selected.species]-required.data[[2]][,selected.species])^2+(out1[,selected.species]-required.data[[3]][,selected.species])^2+(out1[,selected.species]-required.data[[4]][,selected.species])^2)/3 			#22.09.2014
	error.matrix<-(2*alpha.var+1/2)*log(1+sqrt.dif/(2*beta.var))			#22.09.2014
	
	error.matrix[error.matrix=="NaN"]<-Inf
	error.matrix<-sum(error.matrix)

	exponential<--1*sum(error.matrix)

	if(exponential==-Inf)exponential<--1e+300
	}
return(exponential)
}

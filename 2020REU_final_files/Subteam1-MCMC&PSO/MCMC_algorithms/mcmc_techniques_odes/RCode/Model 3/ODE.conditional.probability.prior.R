ode.conditional.probability<-function(ks,required.data,alpha.var,beta.var)
{
library(deSolve)

initial.data<-required.data[[1]]
data.matrix1<-required.data[[2]]
data.matrix2<-required.data[[3]]
data.matrix3<-required.data[[4]]
sigma.data<-required.data[[5]]

parms<-c(KCAT=ks[1],K1=ks[2],k5nc=ks[3],k5cn=ks[4],k4nc=ks[5],k4cn=ks[6],k2a=ks[7],k2d=ks[8],k3=ks[9],k6d=ks[10],k6a=ks[11],Vmax7=ks[12],K7=ks[13])
times<-c(0,0.25,0.5,1,1.5,3,4.5,6,7.5,9,10.5,12,18,24,36,48,100,250,500)
xstart<-c(initial.data['receptors'],initial.data['R_Smad_cyt'],initial.data['R_Smad_P_cyt'],initial.data['Smad4_cyt'],initial.data['R_Smad_P_Smad4_cyt'],initial.data['R_Smad_P_Smad4_nuc'],initial.data['R_Smad_nuc'],initial.data['R_Smad_P_nuc'],initial.data['Smad4_nuc'],initial.data['Pi'])


out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname =  "smad", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes"),error=function(e) -Inf)
out1<-data.frame(out1)
cols<-ncol(out1)
sim.length<-nrow(out1)
tim.length<-length(times)


if(sim.length!=tim.length){
	exponential<--1e+300
	#exponential<--Inf
}else if(!all(times %in% out1[,1])){
	exponential<--1e+300
	#exponential<--Inf
}else{
	#print('computing the loglikelihood')
	selected.species<-c('receptors','R_Smad_cyt','R_Smad_P_cyt','Smad4_cyt','R_Smad_P_Smad4_cyt','R_Smad_P_Smad4_nuc','R_Smad_nuc','R_Smad_P_nuc','Smad4_nuc','Pi')
	#------------------------------ Computation with data points ----------------------------------------
	sqrt.dif<-((out1[,selected.species]-required.data[[2]][,selected.species])^2+(out1[,selected.species]-required.data[[3]][,selected.species])^2+(out1[,selected.species]-required.data[[4]][,selected.species])^2)/3 			#22.09.2014
	error.matrix<-(2*alpha.var+1/2)*log(1+sqrt.dif/(2*beta.var))			#22.09.2014
	
	#error.matrix[error.matrix=="NaN"]<-Inf
	error.matrix[error.matrix=="NaN"]<-1e300
	error.matrix<-sum(error.matrix)

	exponential<--1*sum(error.matrix)
}


return(exponential)
}

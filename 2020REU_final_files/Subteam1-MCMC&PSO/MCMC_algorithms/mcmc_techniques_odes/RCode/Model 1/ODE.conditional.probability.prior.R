ode.conditional.probability<-function(ks,required.data,alpha.var,beta.var)
{
library(deSolve)

initial.data<-required.data[[1]]
data.matrix1<-required.data[[2]]
data.matrix2<-required.data[[3]]
data.matrix3<-required.data[[4]]
sigma.data<-required.data[[5]]

parms<-c(V1=ks[1],n=ks[2],Ki=ks[3],K1=ks[4],V2=ks[5],K2=ks[6],k3=ks[7],K3=ks[8],k4=ks[9],K4=ks[10],V5=ks[11],K5=ks[12],V6=ks[13],K6=ks[14],k7=ks[15],K7=ks[16],k8=ks[17],K8=ks[18],V9=ks[19],K9=ks[20],V10=ks[21],K10=ks[22])

#print(parms)

times<-c(0,0.25,0.5,1,1.5,3,4.5,6,7.5,9,10.5,12,18,24,36,48,100)

xstart<-c(initial.data["MKKK"],initial.data["MKKK_P"],initial.data["MKK"],initial.data["MKK_P"],initial.data["MKK_PP"],initial.data["MAPK"],initial.data["MAPK_P"],initial.data["MAPK_PP"])

############## original solver######################3
###### function for problems with ode solver
out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk_feedback", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes"),error=function(e) -Inf)

out1<-data.frame(out1)
cols<-ncol(out1)

sim.length<-nrow(out1)
tim.length<-length(times)

if(sim.length!=tim.length){
	exponential<--1e+300
	#exponential<--Inf
}else if(!all(times%in% out1[,1] )){
	exponential<--1e+300
#	exponential<--Inf
}else{
	selected.species<-c("MKKK","MKKK_P","MKK","MKK_P","MKK_PP","MAPK","MAPK_P","MAPK_PP")

	sqrt.dif<-((out1[,selected.species]-required.data[[2]][,selected.species])^2+(out1[,selected.species]-required.data[[3]][,selected.species])^2+(out1[,selected.species]-required.data[[4]][,selected.species])^2)/3 			#22.09.2014
	error.matrix<-(2*alpha.var+1/2)*log(1+sqrt.dif/(2*beta.var))			#22.09.2014

	error.matrix[error.matrix=="NaN"]<-1e300
	#error.matrix[error.matrix=="NaN"]<-Inf
	error.matrix<-sum(error.matrix)
	exponential<--1*sum(error.matrix)

	}


return(exponential)
}

#setwd("/R codes model/Model 3/0.5")
#source('solver.all.points.R')

setwd("/R codes model/Model 3")
system("R CMD SHLIB smad.c") #ozi
dyn.load("smad.so") #ozi
## solution B=Ax
library('limSolve')
library(deSolve)
setwd("/R codes model/Model 3/0.5")

proteins<-read.csv("noise data 0.5%.csv", row.names=1)

par.k<-read.csv("set parameter solutions_egf path steady state daniel.csv", row.names=1)
par.k<-log(par.k[,1])

point.solver<-function(ks)
{

	parms<-c(KCAT=ks[1],K1=ks[2],k5nc=ks[3],k5cn=ks[4],k4nc=ks[5],k4cn=ks[6],k2a=ks[7],k2d=ks[8],k3=ks[9],k6d=ks[10],k6a=ks[11],Vmax7=ks[12],K7=ks[13])
	times<-c(0,0.25,0.5,1,1.5,3,4.5,6,7.5,9,10.5,12,18,24,36,48,100)

	final.out<-c()
	final.proteins<-c()

	for(i in 1:(length(times)-1))
	{
		print(i)

		xstart<-c(proteins[proteins[,'time']==times[i],'receptors'],proteins[proteins[,'time']==times[i],'R_Smad_cyt'],proteins[proteins[,'time']==times[i],'R_Smad_P_cyt'],proteins[proteins[,'time']==times[i],'Smad4_cyt'],proteins[proteins[,'time']==times[i],'R_Smad_P_Smad4_cyt'],proteins[proteins[,'time']==times[i],'R_Smad_P_Smad4_nuc'],proteins[proteins[,'time']==times[i],'R_Smad_nuc'],proteins[proteins[,'time']==times[i],'R_Smad_P_nuc'],proteins[proteins[,'time']==times[i],'Smad4_nuc'],proteins[proteins[,'time']==times[i],'Pi'])
		out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "smad", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes"),error=function(e) Inf)
		colnames(out1)<-c('time','receptors','R_Smad_cyt','R_Smad_P_cyt','Smad4_cyt','R_Smad_P_Smad4_cyt','R_Smad_P_Smad4_nuc','R_Smad_nuc','R_Smad_P_nuc','Smad4_nuc','Pi','Sum')
		final.out<-rbind(final.out,out1[2,c('receptors','R_Smad_cyt','R_Smad_P_cyt','Smad4_cyt','R_Smad_P_Smad4_cyt','R_Smad_P_Smad4_nuc','R_Smad_nuc','R_Smad_P_nuc','Smad4_nuc','Pi')])
		final.proteins<-rbind(final.proteins,proteins[proteins[,'time']==times[i+1],c('receptors','R_Smad_cyt','R_Smad_P_cyt','Smad4_cyt','R_Smad_P_Smad4_cyt','R_Smad_P_Smad4_nuc','R_Smad_nuc','R_Smad_P_nuc','Smad4_nuc','Pi')])
	}

	print(final.out)
	print(is.numeric(xstart))
	final.out<-final.out[,2:(ncol(final.out)-1)] ### deleting time and sum column
	dif.out<-final.out-final.proteins
	chi.square<-sum(dif.out^2)
	return(chi.square)
}

k.counter<-0

min.function<-function(x)
{
	#k0=0.433  #just here in mcmc i dont use this parameter
	#k.counter<-k.counter+1
	print('hola')
	point.solver(ks=exp(x))	
	#ks=exp(x)
}




ui.m<-diag(1,13)
ci.v<-rep(1e-7,13)
P<-constrOptim(rep(1,13),min.function,NULL,ui=ui.m,ci=ci.v)

parameters<-P$par
names(parameters)<-c('KCAT','K1','k5nc','k5cn','k4nc','k4cn','k2a','k2d','k3','k6d','k6a','Vmax7','K7')

write.csv(data.matrix(parameters), file="set parameter solutions_egf path_SolverAllPoints.csv")

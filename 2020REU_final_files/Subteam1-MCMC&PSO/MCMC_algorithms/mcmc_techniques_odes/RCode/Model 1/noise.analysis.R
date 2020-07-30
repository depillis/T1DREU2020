#setwd("/R codes model1/Model 1")
#source("noise.analysis.R")
system("R CMD SHLIB mapkk.c") #ozi
dyn.load("mapkk.so") #ozi
library(deSolve)

### initial data
xstart<-c(E1=3e-5,KKK=3e-3,E1_KKK=0,E2=3e-4,P_KKK=0,E2_P_KKK=0,KK=1.2,P_KKK_KK=0,P_KK=0,KKPase=3e-4,KKPase_P_KK=0,P_KKK_P_KK=0,PP_KK=0,KKPase_PP_KK=0,K=1.2,PP_KK_K=0,P_K=0,KPase=0.12,KPase_P_K=0,PP_KK_P_K=0,PP_K=0,KPase_PP_K=0)

##compartment=4e-12

parms<-c(a1=1000,d1=150,k1=150,a2=1000,d2=150,k2=150,a3=1000,d3=150,k3=150,a4=1000,d4=150,k4=150,a5=1000,d5=150,k5=150,a6=1000,d6=150,k6=150,a7=1000,d7=150,k7=150,a8=1000,d8=150,k8=150,a9=1000,d9=150,k9=150,a10=1000,d10=150,k10=150,compartment=4e-12)
times<-seq(0, 48, by=0.01 )

out1<-ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes") ##ozi

pdf('possible system.pdf')
protein.name<-colnames(out1)
for(i in 2:(ncol(out1)-1))
{
	plot(out1[,1],out1[,i],ylab=protein.name[i], xlab="Time",type='line')

}

dev.off()

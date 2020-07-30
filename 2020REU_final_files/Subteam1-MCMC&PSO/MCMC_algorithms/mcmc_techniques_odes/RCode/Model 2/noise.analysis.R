#setwd("/20130225_noise analysis2")
#source("noise.analysis.R")
system("R CMD SHLIB jak.c") #ozi
dyn.load("jak.so") #ozi
library(deSolve)

### initial data
xstart<-c(JAK2=2.0,pJAK2=0.0,Epo=50,EpoR=1.0,pEpoR=0.0,SHP1=10.7991,mSHP1=0.0,Delay01_mSHP1=0.0,Delay02_mSHP1=0.0,Delay03_mSHP1=0.0,Delay04_mSHP1=0.0,Delay05_mSHP1=0.0,Delay06_mSHP1=0.0,Delay07_mSHP1=0.0,Delay08_mSHP1=0.0,actSHP1=0.0,SOS=2.5101,mSOS=0.0,Raf=3.7719,pRaf=0.0,MEK2=11.0,pMEK2=0.0,MEK1=24.0,pMEK1=0.0,ppMEK2=0.0,ppMEK1=0.0,ERK1=7.0,pERK1=0.0,ERK2=21.0,pERK2=0.0,ppERK1=0.0,ppERK2=0.0,pSOS=0.0)
##compartment=4e-12

parms<-c(k1=0.0122149,k2=3.15714,k3=0.408408,k4=0.408408,k13=0.0248773,k14=1.19995,k15=0.368384,k16=0.10271,k17=15.5956,k18=0.144515,k19=0.374228,k20=3.11919,k21=0.687193,k22=215.158,k23=667.957,k24=0.130937,k26=0.0732724,k28=2.4927,k29=2.44361,k32=59.5251,k33=53.0816,k36=39.0886,k38=3.00453,k40=5122.68,k42=0.124944,cell=1)
times<-seq(0, 60, by=0.01 )

out1<-ode(xstart,times,func = "derivs",parms = parms, dllname = "jak", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes") ##ozi

pdf('possible system.pdf')
protein.name<-colnames(out1)
for(i in 2:(ncol(out1)-1))
{
	plot(out1[,1],out1[,i],ylab=protein.name[i], xlab="Time",type='line')

}

dev.off()



times<-c(0,0.25,0.5,1,1.5,3,4.5,6,7.5,9,10.5,12,18,24,36,48,100)
out1<-ode(xstart,times,func = "derivs",parms = parms, dllname = "jak", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes") ##ozi

pdf('possible system_less time point.pdf')
protein.name<-colnames(out1)
for(i in 2:(ncol(out1)-1))
{
	plot(out1[,1],out1[,i],ylab=protein.name[i], xlab="Time[min]",type='line',xlim=c(0,100))
}
dev.off()


write.csv(out1[,1:(ncol(out1)-1)],file='data without noise.csv')



sim.data<-out1[,1:(ncol(out1)-1)]
new.sim.data<-sim.data
new.sim.data2<-sim.data
new.sim.data3<-sim.data
cols<-ncol(sim.data)
rows<-nrow(sim.data)

sd.vector<-c(0.001,0.005,0.01,0.05,0.1,0.15,0.25)

for(j in 1:length(sd.vector))
{
	for(i in 2:cols)
	{
		amp.data<-max(sim.data[,i])-min(sim.data[,i])
		noise.data<-rnorm(rows,sd=sd.vector[j])
		new.sim.data[,i]<-sim.data[,i]+amp.data*noise.data	

		noise.data<-rnorm(rows,sd=sd.vector[j])
		new.sim.data2[,i]<-sim.data[,i]+amp.data*noise.data

		noise.data<-rnorm(rows,sd=sd.vector[j])
		new.sim.data2[,i]<-sim.data[,i]+amp.data*noise.data
		
	}
	write.csv(new.sim.data,paste('noise data ',(100*sd.vector[j]),'%_r1.csv',sep=''),sep="\t")
	write.csv(new.sim.data2,paste('noise data ',(100*sd.vector[j]),'%_r2.csv',sep=''),sep="\t")
	write.csv(new.sim.data3,paste('noise data ',(100*sd.vector[j]),'%_r3.csv',sep=''),sep="\t")
}

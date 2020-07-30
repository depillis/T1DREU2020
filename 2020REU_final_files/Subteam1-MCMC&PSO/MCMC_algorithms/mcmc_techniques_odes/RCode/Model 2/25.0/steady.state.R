#setwd("/R codes model/Model 2/25.0")
#source('steady.state.R')

## solution B=Ax
library('limSolve')

dproteins<-read.csv("derivation data 16 time course.smooth.csv",row.names=1)
proteins<-read.csv("noise data 25%.csv", row.names=1)

JAK2<-proteins[proteins[,"time"]==0,"JAK2"]
dJAK2<-0

pJAK2<-proteins[proteins[,"time"]==0,"pJAK2"]
dpJAK2<-0

Epo<-proteins[proteins[,"time"]==0,"Epo"]
dEpo<-0

EpoR<-proteins[proteins[,"time"]==0,"EpoR"]
dEpoR<-0

pEpoR<-proteins[proteins[,"time"]==0,"pEpoR"]
dpEpoR<-0

SHP1<-proteins[proteins[,"time"]==0,"SHP1"]
dSHP1<-0

mSHP1<-proteins[proteins[,"time"]==0,"mSHP1"]
dmSHP1<-0

Delay01_mSHP1<-proteins[proteins[,"time"]==0,"Delay01_mSHP1"]
dDelay01_mSHP1<-0

Delay02_mSHP1<-proteins[proteins[,"time"]==0,"Delay02_mSHP1"]
dDelay02_mSHP1<-0

Delay03_mSHP1<-proteins[proteins[,"time"]==0,"Delay03_mSHP1"]
dDelay03_mSHP1<-0

Delay04_mSHP1<-proteins[proteins[,"time"]==0,"Delay04_mSHP1"]
dDelay04_mSHP1<-0

Delay05_mSHP1<-proteins[proteins[,"time"]==0,"Delay05_mSHP1"]
dDelay05_mSHP1<-0

Delay06_mSHP1<-proteins[proteins[,"time"]==0,"Delay06_mSHP1"]
dDelay06_mSHP1<-0

Delay07_mSHP1<-proteins[proteins[,"time"]==0,"Delay07_mSHP1"]
dDelay07_mSHP1<-0

Delay08_mSHP1<-proteins[proteins[,"time"]==0,"Delay08_mSHP1"]
dDelay08_mSHP1<-0

actSHP1<-proteins[proteins[,"time"]==0,"actSHP1"]
dactSHP1<-0

SOS<-proteins[proteins[,"time"]==0,"SOS"]
dSOS<-0

mSOS<-proteins[proteins[,"time"]==0,"mSOS"]
dmSOS<-0

Raf<-proteins[proteins[,"time"]==0,"Raf"]
dRaf<-0

pRaf<-proteins[proteins[,"time"]==0,"pRaf"]
dpRaf<-0

MEK2<-proteins[proteins[,"time"]==0,"MEK2"]
dMEK2<-0

pMEK2<-proteins[proteins[,"time"]==0,"pMEK2"]
dpMEK2<-0

MEK1<-proteins[proteins[,"time"]==0,"MEK1"]
dMEK1<-0

pMEK1<-proteins[proteins[,"time"]==0,"pMEK1"]
dpMEK1<-0

ppMEK2<-proteins[proteins[,"time"]==0,"ppMEK2"]
dppMEK2<-0

ppMEK1<-proteins[proteins[,"time"]==0,"ppMEK1"]
dppMEK1<-0

ERK1<-proteins[proteins[,"time"]==0,"ERK1"]
dERK1<-0

pERK1<-proteins[proteins[,"time"]==0,"pERK1"]
dpERK1<-0

ERK2<-proteins[proteins[,"time"]==0,"ERK2"]
dERK2<-0

pERK2<-proteins[proteins[,"time"]==0,"pERK2"]
dpERK2<-0

ppERK1<-proteins[proteins[,"time"]==0,"ppERK1"]
dppERK1<-0

ppERK2<-proteins[proteins[,"time"]==0,"ppERK2"]
dppERK2<-0

pSOS<-proteins[proteins[,"time"]==0,"pSOS"]
dpSOS<-0


min.function<-function(x)
{
	k1=x[1]    	##  JAK2_phosphorylation_by_Epo   		##
	k2=x[2]    	##  EpoR_phosphorylation_by_pJAK2 		##
	k3=x[3]    	##  SHP1_activation_by_pEpoR      		##
	k4=x[4]    	##  SHP1_delay			 		##
	k13=x[5]  	##  actSHP1_deactivation			##
	k14=x[6]  	##  pEpoR_desphosphorylation_by_actSHP1 	##
	k15=x[7] 	##  pJAK2_dephosphorylation_by_actSHP		##
	k16=x[8]  	##  SOS_recruitment_by_pEpoR			##
	k17=x[9]  	##  mSOS_release_from_membrane			##
	k18=x[10] 	##  mSOS_induced_Raf_phosphorylation		##
	k19=x[11] 	##  pRaf_dephosphorylation			##
	k20=x[12]	##  First_MEK2_phosphorylation_by_pRaf		##
	k21=x[13]	##  First_MEK1_phosphorylation_by_pRaf		##
	k22=x[14]	##  Second_MEK2_phosphorylation_by_pRaf		##
	k23=x[15]	##  Second_MEK1_phosphorylation_by_pRaf		##
	k24=x[16]	##  First_MEK_desphosphorylation		##
	k26=x[17]	##  Second_MEK_desphosphorylation		##
	k28=x[18]	##  First_ERK1_phosphorylation_by_ppMEK		##
	k29=x[19]	##  First_ERK2_phosphorylation_by_ppMEK		##
	k32=x[20]	##  Second_ERK1_phosphorylation_by_ppMEK	##
	k33=x[21]	##  Second_ERK2_phosphorylation_by_ppMEK	##
	k36=x[22]	##  First_ERK_dephosphorylation			##
	k38=x[23]	##  Second_ERK_dephosphorylation		##
	k40=x[24]	##  ppERK_neg_feedback_on_mSOS			##
	k42=x[25]	##  pSOS_dephosphorylation			##
	cell=1

	sum((dJAK2-(-k1*JAK2*Epo*cell+k15*pJAK2*actSHP1*cell))^2+(dpJAK2-(k1*JAK2*Epo*cell-k15*pJAK2*actSHP1*cell))^2+(dEpo-(0))^2+(dEpoR-(-k2*EpoR*pJAK2*cell+k14*pEpoR*actSHP1*cell))^2+(dpEpoR-(k2*EpoR*pJAK2*cell-k14*pEpoR*actSHP1*cell))^2+
	(dSHP1-(-k3*SHP1*pEpoR*cell+k13*actSHP1*cell))^2+(dmSHP1-(k3*SHP1*pEpoR*cell-k4*mSHP1*cell))^2+(dDelay01_mSHP1-(k4*mSHP1*cell-k4*Delay01_mSHP1*cell))^2+(dDelay02_mSHP1-(k4*Delay01_mSHP1*cell-k4*Delay02_mSHP1*cell))^2+(dDelay03_mSHP1-(k4*Delay02_mSHP1*cell-k4*Delay03_mSHP1*cell))^2+(dDelay04_mSHP1-(k4*Delay03_mSHP1*cell-k4*Delay04_mSHP1*cell))^2+(dDelay05_mSHP1-(k4*Delay04_mSHP1*cell-k4*Delay05_mSHP1*cell))^2+(dDelay06_mSHP1-(k4*Delay05_mSHP1*cell-k4*Delay06_mSHP1*cell))^2+(dDelay07_mSHP1-(k4*Delay06_mSHP1*cell-k4*Delay07_mSHP1*cell))^2+(dDelay08_mSHP1-(k4*Delay07_mSHP1*cell-k4*Delay08_mSHP1*cell))^2+(dactSHP1-(k4*Delay08_mSHP1*cell-k13*actSHP1*cell))^2+(dSOS-(-k16*SOS*pEpoR*cell+k17*mSOS*cell+k42*pSOS*cell))^2+(dmSOS-(k16*SOS*pEpoR*cell-k17*mSOS*cell-k40*mSOS*ppERK1*cell-k40*mSOS*ppERK2*cell))^2+(dRaf-(-k18*Raf*mSOS*cell+k19*pRaf*cell))^2+(dpRaf-(k18*Raf*mSOS*cell-k19*pRaf*cell))^2+(dMEK2-(-k20*MEK2*pRaf*cell+k26*pMEK2*cell))^2+(dpMEK2-(k20*MEK2*pRaf*cell-k22*pMEK2*pRaf*cell+k24*ppMEK2*cell-k26*pMEK2*cell))^2+(dMEK1-(-k21*MEK1*pRaf*cell+k26*pMEK1*cell))^2+(dpMEK1-(k21*MEK1*pRaf*cell-k23*pMEK1*pRaf*cell+k24*ppMEK1*cell-k26*pMEK1*cell))^2+(dppMEK2-(k22*pMEK2*pRaf*cell-k24*ppMEK2*cell))^2+(dppMEK1-(k23*pMEK1*pRaf*cell-k24*ppMEK1*cell))^2+(dERK1-(-k28*ERK1*ppMEK2*cell-k28*ERK1*ppMEK1*cell+k38*pERK1*cell))^2+(dpERK1-(k28*ERK1*ppMEK2*cell+k28*ERK1*ppMEK1*cell-k32*pERK1*ppMEK2*cell-k32*pERK1*ppMEK1*cell+k36*ppERK1*cell-k38*pERK1*cell))^2+(dERK2-(-k29*ERK2*ppMEK2*cell-k29*ERK2*ppMEK1*cell+k38*pERK2*cell))^2+(dpERK2-(k29*ERK2*ppMEK2*cell+k29*ERK2*ppMEK1*cell-k33*pERK2*ppMEK2*cell-k33*pERK2*ppMEK1*cell+k36*ppERK2*cell-k38*pERK2*cell))^2+(dppERK1-(k32*pERK1*ppMEK2*cell+k32*pERK1*ppMEK1*cell-k36*ppERK1*cell))^2+(dppERK2-(k33*pERK2*ppMEK2*cell+k33*pERK2*ppMEK1*cell-k36*ppERK2*cell))^2+(dpSOS-(k40*mSOS*ppERK1*cell+k40*mSOS*ppERK2*cell-k42*pSOS*cell))^2)
}


ui.m<-diag(1,25)
ci.v<-rep(1e-7,25)
P<-constrOptim(rep(1,25),min.function,NULL,ui=ui.m,ci=ci.v)

parameters<-P$par

names(parameters)<-c('k1','k2','k3','k4','k13','k14','k15','k16','k17','k18','k19','k20','k21','k22','k23','k24','k26','k28','k29','k32','k33','k36','k38','k40','k42')
write.csv(parameters,file="set parameter solutions_egf path steady state daniel.csv")

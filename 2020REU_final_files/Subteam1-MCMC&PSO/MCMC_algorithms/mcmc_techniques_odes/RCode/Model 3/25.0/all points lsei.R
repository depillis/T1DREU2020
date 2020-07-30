#setwd("/R codes model/Model 3/25.0")
#source('all points lsei.R')

## solution B=Ax
library('limSolve')


dproteins<-read.csv("derivation data 16 time course.smooth.csv", row.names=1)


proteins<-read.csv("noise data 25%.csv", row.names=1)

receptors<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"receptors"]
receptors<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"receptors"]

R_Smad_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_cyt"]
R_Smad_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_cyt"]

R_Smad_P_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_cyt"]
R_Smad_P_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_cyt"]

Smad4_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"Smad4_cyt"]
Smad4_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"Smad4_cyt"]

R_Smad_P_Smad4_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_Smad4_cyt"]
R_Smad_P_Smad4_cyt<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_Smad4_cyt"]

R_Smad_P_Smad4_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_Smad4_nuc"]
R_Smad_P_Smad4_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_Smad4_nuc"]

R_Smad_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_nuc"]
R_Smad_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_nuc"]

R_Smad_P_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_nuc"]
R_Smad_P_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"R_Smad_P_nuc"]

Smad4_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"Smad4_nuc"]
Smad4_nuc<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"Smad4_nuc"]

Pi<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"Pi"]
Pi<-proteins[proteins[,"time"]>=0&proteins[,"time"]<500,"Pi"]

min.function<-function(x)
{
	KCAT=x[1]
	K1=x[2]
	k5nc=x[3]
	k5cn=x[4]
	k4nc=x[5]
	k4cn=x[6]
	k2a=x[7]
	k2d=x[8]
	k3=x[9]
	k6d=x[10]
	k6a=x[11]
	Vmax7=x[12]
	K7=x[13]

	sum(+(dreceptors-(-1*receptors/90))^2+(dR_Smad_cyt-(-KCAT*receptors*R_Smad_cyt/(K1+R_Smad_cyt)+k5nc*R_Smad_nuc-k5cn*R_Smad_cyt))^2+(dR_Smad_P_cyt-(KCAT*receptors*R_Smad_cyt/(K1+R_Smad_cyt)-k2a*R_Smad_P_cyt*Smad4_cyt+k2d*R_Smad_P_Smad4_cyt))^2+(dSmad4_cyt-(k4nc*Smad4_nuc-k4cn*Smad4_cyt-k2a*R_Smad_P_cyt*Smad4_cyt+k2d*R_Smad_P_Smad4_cyt))^2+(dR_Smad_P_Smad4_cyt-(k2a*R_Smad_P_cyt*Smad4_cyt-k2d*R_Smad_P_Smad4_cyt-k3*R_Smad_P_Smad4_cyt))^2+(dR_Smad_P_Smad4_nuc-(k3*R_Smad_P_Smad4_cyt-k6d*R_Smad_P_Smad4_nuc+k6a*Smad4_nuc*R_Smad_P_nuc))^2+(dR_Smad_nuc-(-k5nc*R_Smad_nuc+k5cn*R_Smad_cyt+Vmax7*R_Smad_P_nuc/(K7+R_Smad_P_nuc)))^2+(dR_Smad_P_nuc-(k6d*R_Smad_P_Smad4_nuc-k6a*Smad4_nuc*R_Smad_P_nuc-Vmax7*R_Smad_P_nuc/(K7+R_Smad_P_nuc)))^2+(dSmad4_nuc-(-k4nc*Smad4_nuc+k4cn*Smad4_cyt+k6d*R_Smad_P_Smad4_nuc-k6a*Smad4_nuc*R_Smad_P_nuc))^2+(dPi-(Vmax7*R_Smad_P_nuc/(K7+R_Smad_P_nuc)))^2)
}


ui.m<-diag(1,13)
ci.v<-rep(1e-7,13)
P<-constrOptim(rep(1,13),min.function,NULL,ui=ui.m,ci=ci.v)

parameters<-P$par
names(parameters)<-c('KCAT','K1','k5nc','k5cn','k4nc','k4cn','k2a','k2d','k3','k6d','k6a','Vmax7','K7')

write.csv(data.matrix(parameters), file="set parameter solutions_egf path all points.csv")

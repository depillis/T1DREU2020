#path<-'/R code model/Model 5/'
setwd(path)
#source('steady.state.R')

## solution B=Ax
setwd(path)
library('limSolve')


#dproteins<-read.csv("derivation data 16 time course.smooth.csv", row.names=1)
proteins<-data.matrix(read.csv("initial.species.values.csv",row.names=1))


Rec <- proteins["Rec",1]
##Rec_i <- proteins["Rec_i",1]  # this is a parameter, it does not have intial value
IL13_Rec<- proteins["IL13_Rec",1]
p_IL13_Rec <- proteins["p_IL13_Rec",1]
p_IL13_Rec_i <- proteins["p_IL13_Rec_i",1]
JAK2 <- proteins["JAK2",1]
pJAK2 <- proteins["pJAK2",1]
STAT5 <- proteins["STAT5",1]
pSTAT5 <- proteins["pSTAT5",1]
SOCS3mRNA <- proteins["SOCS3mRNA",1]
DecoyR <- proteins["DecoyR",1]
IL13_DecoyR <- proteins["IL13_DecoyR",1]
SOCS3 <- proteins["SOCS3",1]
CD274mRNA <- proteins["CD274mRNA",1]


dRec <- 0
dRec_i <- 0
dIL13_Rec <- 0
dp_IL13_Rec <- 0
dp_IL13_Rec_i <- 0
dJAK2 <- 0
dpJAK2 <- 0
dSTAT5 <- 0
dpSTAT5 <- 0
dSOCS3mRNA <- 0
dDecoyR <- 0
dIL13_DecoyR <- 0
dSOCS3 <- 0
dCD274mRNA <- 0


observable.data<-data.matrix(read.csv(paste(path,'/Observable.data.csv',sep=''),sep=';'))

#sigmaY1TimR = sd(observable.data[observable.data[,'Model']==1,'IR1_P'])
#sigmaY2Step = sd(observable.data[observable.data[,'Model']==1,'IRS1_P'])
#sigmaY2TimR = sd(observable.data[observable.data[,'Model']==2,'IRS1_P'])
#sigmaYDosR = sd(observable.data[observable.data[,'Model']!=2&observable.data[,'Model']!=1,'IRS1_P'])


min.function<-function(x)
{
#kinetics parameters

  CD274mRNA_production = x[1]
  DecoyR_binding = x[2]
  JAK2_p_inhibition = x[3]
  JAK2_phosphorylation = x[4]
  Kon_IL13Rec = x[5]
  Rec_intern = x[6]
  Rec_phosphorylation = x[7]
  Rec_recycle = x[8]
  SOCS3_accumulation = x[9]
  SOCS3_degradation = x[10]
  SOCS3_translation = x[11]
  SOCS3mRNA_production = x[12]
  STAT5_phosphorylation = x[13]
  pJAK2_dephosphorylation = x[14]
  pRec_degradation = x[15]
  pRec_intern = x[16]
  pSTAT5_dephosphorylation = x[17]
#initial species value as parameter
  Rec_i = x[18] 
#scalar parameters
  scale_CD274mRNA_obs = x[19]
  scale_IL13_cell_obs = x[20]
  scale_SOCS3mRNA_obs = x[21]
  scale_pIL4Ra_obs = x[22]
  scale_pJAK2_obs = x[23]

#k_IRP_1Step = x[15]
#k_IRSiP_1Step = x[16]
#k_IRSiP_2Step = x[17]
#k_IRSiP_DosR = x[18]

SHP1 <- 91
il13_level <- 0

IL13_cell <- IL13_Rec + IL13_DecoyR + p_IL13_Rec + p_IL13_Rec_i
pIL4Ra <- p_IL13_Rec + p_IL13_Rec_i


#sum((dMKKK-(V2*MKKK_P/(K2+MKKK_P)-V1*MKKK/((1+(MAPK_PP/Ki)^n)*(K1+MKKK))))^2+(dMKKK_P-(V1*MKKK/((1+(MAPK_PP/Ki)^n)*(K1+MKKK)) - V2*MKKK_P/(K2+MKKK_P)))^2+(dMKK-(V6*MKK_P/(K6+MKK_P) - k3*MKKK_P*MKK/(K3+MKK)))^2+(dMKK_P-(k3*MKKK_P*MKK/(K3+MKK) + V5*MKK_PP/(K5+MKK_PP) - k4*MKKK_P*MKK_P/(K4+MKK_P)  - V6*MKK_P/(K6+MKK_P)))^2+(dMKK_PP-(k4*MKKK_P*MKK_P/(K4+MKK_P)  - V5*MKK_PP/(K5+MKK_PP)))^2+(dMAPK-(V10*MAPK_P/(K10+MAPK_P)  - k7*MKK_PP*MAPK/(K7+MAPK)))^2+(dMAPK_P-(k7*MKK_PP*MAPK/(K7+MAPK) + V9*MAPK_PP/(K9+MAPK_PP) - k8*MKK_PP*MAPK_P/(K8+MAPK_P) - V10*MAPK_P/(K10+MAPK_P)))^2+(dMAPK_PP-(k8*MKK_PP*MAPK_P/(K8+MAPK_P)  - V9*MAPK_PP/(K9+MAPK_PP)))^2)

return(sum( ( dRec - (Rec_i*Rec_recycle - Rec*Rec_intern - il13_level*Kon_IL13Rec*Rec) )^2
	+( dRec_i - (Rec*Rec_intern - Rec_i*Rec_recycle) )^2
	+( dIL13_Rec - (il13_level*Kon_IL13Rec*Rec - IL13_Rec*Rec_phosphorylation*pJAK2) )^2
 	+( dp_IL13_Rec - (IL13_Rec*Rec_phosphorylation*pJAK2 - pRec_intern*p_IL13_Rec) )^2
	+( dp_IL13_Rec_i - (pRec_intern*p_IL13_Rec - pRec_degradation*p_IL13_Rec_i) )^2
	+( dJAK2 - (SHP1*pJAK2*pJAK2_dephosphorylation - (JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)) )^2
	+( dpJAK2 - ((JAK2*JAK2_phosphorylation*p_IL13_Rec)/(JAK2_p_inhibition*SOCS3 + 1) - SHP1*pJAK2*pJAK2_dephosphorylation + (IL13_Rec*JAK2*JAK2_phosphorylation)/(JAK2_p_inhibition*SOCS3 + 1)) )^2
	+( dSTAT5 - (SHP1*pSTAT5*pSTAT5_dephosphorylation - STAT5*STAT5_phosphorylation*pJAK2) )^2
	+( dpSTAT5 - (STAT5*STAT5_phosphorylation*pJAK2 - SHP1*pSTAT5*pSTAT5_dephosphorylation) )^2
	+( dSOCS3mRNA - (SOCS3mRNA_production*pSTAT5) )^2
	+( dDecoyR - (-il13_level*DecoyR*DecoyR_binding) )^2
	+( dIL13_DecoyR - (il13_level*DecoyR*DecoyR_binding) )^2
	+( dSOCS3 - ((SOCS3mRNA*SOCS3_translation)/(SOCS3mRNA + SOCS3_accumulation) - SOCS3*SOCS3_degradation) )^2
	+( dCD274mRNA - (CD274mRNA_production*pSTAT5) )^2
	+ (Rec + IL13_Rec + p_IL13_Rec - observable.data[observable.data[,1]==1&observable.data[,2]==0,'RecSurf_obs'] )^2
	+ (scale_IL13_cell_obs * IL13_cell - observable.data[observable.data[,1]==1&observable.data[,2]==0,'IL13_cell_obs'] )^2
	+ (scale_pIL4Ra_obs * pIL4Ra - observable.data[observable.data[,1]==1&observable.data[,2]==0,'pIL4Ra_obs'] )^2
	+ (scale_pJAK2_obs * pJAK2 - observable.data[observable.data[,1]==1&observable.data[,2]==0,'pJAK2_obs'] )^2
	+ (scale_SOCS3mRNA_obs * SOCS3mRNA - observable.data[observable.data[,1]==1&observable.data[,2]==0,'SOCS3mRNA_obs'] )^2
	+ (scale_CD274mRNA_obs * CD274mRNA - observable.data[observable.data[,1]==1&observable.data[,2]==0,'CD274mRNA_obs'] )^2
	+ (SOCS3 - observable.data[observable.data[,1]==1&observable.data[,2]==0,'SOCS3_obs'] )^2
	+ (pSTAT5 - observable.data[observable.data[,1]==1&observable.data[,2]==0,'pSTAT5_obs'] )^2 ) )
}


ui.m<-diag(1,23)
ci.v<-rep(1e-7,23)
P<-constrOptim(rep(1,23),min.function,NULL,ui=ui.m,ci=ci.v)

#parameters<-c(P$par,sigmaY1TimR,sigmaY2Step,sigmaY2TimR,sigmaYDosR)
parameters<-c(P$par)


#names(parameters)<-c('V1','n','Ki','K1','V2','K2','k3','K3','k4','K4','V5','K5','V6','K6','k7','K7','k8','K8','V9','K9','V10','K10')
#names(parameters)<-c('k1a','k1aBasic','k1b','k1c','k1d','k1e','k1f','k1g','k1r','k21','k22','k3','km2','km3','k_IRP_1Step','k_IRSiP_1Step','k_IRSiP_2Step','k_IRSiP_DosR','sigmaY1TimR','sigmaY2Step','sigmaY2TimR','sigmaYDosR' )
names(parameters)<-c("CD274mRNA_production", "DecoyR_binding",
 "JAK2_p_inhibition", "JAK2_phosphorylation",
  "Kon_IL13Rec", "Rec_intern", "Rec_phosphorylation", 
  "Rec_recycle", "SOCS3_accumulation", "SOCS3_degradation", 
  "SOCS3_translation", "SOCS3mRNA_production", 
  "STAT5_phosphorylation", "pJAK2_dephosphorylation", 
  "pRec_degradation", "pRec_intern", "pSTAT5_dephosphorylation", 
  "Rec_i", "scale_CD274mRNA_obs", "scale_IL13_cell_obs", 
  "scale_SOCS3mRNA_obs", "scale_pIL4Ra_obs", "scale_pJAK2_obs")


#write.csv(data.matrix(parameters),file="set parameter solutions_egf path steady state.daniel.csv")
write.csv(parameters,file=paste(path,"realistic/set parameter solutions_egf path steady state daniel.csv",sep=''))



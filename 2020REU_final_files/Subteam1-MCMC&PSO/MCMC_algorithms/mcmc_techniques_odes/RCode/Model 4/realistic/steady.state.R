#path<-'/R codes model/Model 4/'
#source('steady.state.R')

## solution B=Ax
setwd(path)
library('limSolve')


#dproteins<-read.csv("derivation data 16 time course.smooth.csv", row.names=1)
proteins<-data.matrix(read.csv("initial.species.values.csv",row.names=1))


IR <- proteins["IR",1]
IRS <- proteins["IRS",1]
IRSiP <- proteins["IRSiP",1]
IRi <- proteins["IRi",1]
IRiP <- proteins["IRiP",1]
IRins <- proteins["IRins",1]
IRp <- proteins["IRp",1]
XX <- proteins["XX",1]
XXp <- proteins["XXp",1]

dIR <- 0
dIRS <- 0
dIRSiP <- 0
dIRi <- 0
dIRiP <- 0
dIRins <- 0
dIRp <- 0
dXX <- 0
dXXp <- 0


observable.data<-data.matrix(read.csv(paste(path,'/Observable.data.csv',sep=''),sep=';'))

sigmaY1TimR = sd(observable.data[observable.data[,'Model']==1,'IR1_P'])
sigmaY2Step = sd(observable.data[observable.data[,'Model']==1,'IRS1_P'])
sigmaY2TimR = sd(observable.data[observable.data[,'Model']==2,'IRS1_P'])
sigmaYDosR = sd(observable.data[observable.data[,'Model']!=2&observable.data[,'Model']!=1,'IRS1_P'])


min.function<-function(x)
{

k1a = x[1]
k1aBasic = x[2]
k1b = x[3]
k1c = x[4]
k1d = x[5]
k1e = x[6]
k1f = x[7]
k1g = x[8]
k1r = x[9]
k21 = x[10]
k22 = x[11]
k3 = x[12]
km2 = x[13]
km3 = x[14]
k_IRP_1Step = x[15]
k_IRSiP_1Step = x[16]
k_IRSiP_2Step = x[17]
k_IRSiP_DosR = x[18]

insulin <- 0


#sum((dMKKK-(V2*MKKK_P/(K2+MKKK_P)-V1*MKKK/((1+(MAPK_PP/Ki)^n)*(K1+MKKK))))^2+(dMKKK_P-(V1*MKKK/((1+(MAPK_PP/Ki)^n)*(K1+MKKK)) - V2*MKKK_P/(K2+MKKK_P)))^2+(dMKK-(V6*MKK_P/(K6+MKK_P) - k3*MKKK_P*MKK/(K3+MKK)))^2+(dMKK_P-(k3*MKKK_P*MKK/(K3+MKK) + V5*MKK_PP/(K5+MKK_PP) - k4*MKKK_P*MKK_P/(K4+MKK_P)  - V6*MKK_P/(K6+MKK_P)))^2+(dMKK_PP-(k4*MKKK_P*MKK_P/(K4+MKK_P)  - V5*MKK_PP/(K5+MKK_PP)))^2+(dMAPK-(V10*MAPK_P/(K10+MAPK_P)  - k7*MKK_PP*MAPK/(K7+MAPK)))^2+(dMAPK_P-(k7*MKK_PP*MAPK/(K7+MAPK) + V9*MAPK_PP/(K9+MAPK_PP) - k8*MKK_PP*MAPK_P/(K8+MAPK_P) - V10*MAPK_P/(K10+MAPK_P)))^2+(dMAPK_PP-(k8*MKK_PP*MAPK_P/(K8+MAPK_P)  - V9*MAPK_PP/(K9+MAPK_PP)))^2)

sum(( IRins*k1b - IR*k1aBasic + IRp*k1g + IRi*k1r - IR*insulin*k1a - dIR )^2
+ ( IR*k1aBasic - IRins*k1b - IRins*k1c + IR*insulin*k1a - dIRins )^2
+ ( IRins*k1c - IRp*k1d - IRp*k1g - dIRp )^2
+ ( IRp*k1d - IRiP*(k1e + (XXp*k1f)/(XXp + 1)) - dIRiP )^2
+ ( IRiP*(k1e + (XXp*k1f)/(XXp + 1)) - IRi*k1r - dIRi )^2
+ ( IRSiP*km2 - IRS*k21*(IRp + IRiP*k22) - dIRS )^2
+ ( IRS*k21*(IRp + IRiP*k22) - IRSiP*km2 - dIRSiP )^2
+ ( XXp*km3 - IRSiP*XX*k3 - dXX )^2
+ ( IRSiP*XX*k3 - XXp*km3 - dXXp )^2 
+(k_IRP_1Step*(IRiP + IRp)-observable.data[observable.data[,1]==1&observable.data[,2]==0,'IR1_P'])^2
+(IRSiP*k_IRSiP_1Step - observable.data[observable.data[,1]==1&observable.data[,2]==0,'IRS1_P'])^2
## Model 2 ####
+(IRSiP*k_IRSiP_2Step - observable.data[observable.data[,1]==2&observable.data[,2]==0,'IRS1_P'])^2
## Model 3 ####
+(IRSiP*k_IRSiP_DosR - mean(observable.data[observable.data[,1]!=1&observable.data[,1]!=2,'IRS1_P'])  )^2 )
}


ui.m<-diag(1,18)
ci.v<-rep(1e-7,18)
P<-constrOptim(rep(1,18),min.function,NULL,ui=ui.m,ci=ci.v)

parameters<-c(P$par,sigmaY1TimR,sigmaY2Step,sigmaY2TimR,sigmaYDosR)

#names(parameters)<-c('V1','n','Ki','K1','V2','K2','k3','K3','k4','K4','V5','K5','V6','K6','k7','K7','k8','K8','V9','K9','V10','K10')
names(parameters)<-c('k1a','k1aBasic','k1b','k1c','k1d','k1e','k1f','k1g','k1r','k21','k22','k3','km2','km3','k_IRP_1Step','k_IRSiP_1Step','k_IRSiP_2Step','k_IRSiP_DosR','sigmaY1TimR','sigmaY2Step','sigmaY2TimR','sigmaYDosR' )


#write.csv(data.matrix(parameters),file="set parameter solutions_egf path steady state.daniel.csv")
write.csv(parameters,file=paste(path,"realistic/set parameter solutions_egf path steady state daniel.csv",sep=''))

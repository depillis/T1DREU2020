#source('all points lsei.R')

## solution B=Ax
library('limSolve')


dproteins<-read.csv("derivation data 16 time course.smooth.csv", row.names=1)
proteins<-read.csv("noise data 15%.csv", row.names=1)


MKKK<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MKKK"]
dMKKK<-dproteins[dproteins[,"time"]<48,"MKKK"]

MKKK_P<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MKKK_P"]
dMKKK_P<-dproteins[dproteins[,"time"]<48,"MKKK_P"]

MKK<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MKK"]
dMKK<-dproteins[dproteins[,"time"]<48,"MKK"]

MKK_P<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MKK_P"]
dMKK_P<-dproteins[dproteins[,"time"]<48,"MKK_P"]

MKK_PP<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MKK_PP"]
dMKK_PP<-dproteins[dproteins[,"time"]<48,"MKK_PP"]

MAPK<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MAPK"]
dMAPK<-dproteins[dproteins[,"time"]<48,"MAPK"]

MAPK_P<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MAPK_P"]
dMAPK_P<-dproteins[dproteins[,"time"]<48,"MAPK_P"]

MAPK_PP<-proteins[proteins[,"time"]>=0&proteins[,"time"]<48,"MAPK_PP"]
dMAPK_PP<-dproteins[dproteins[,"time"]<48,"MAPK_PP"]





min.function<-function(x)
{
V1=x[1]
n=x[2]
Ki=x[3]
K1=x[4]
V2=x[5]
K2=x[6]
k3=x[7]
K3=x[8]
k4=x[9]
K4=x[10]
V5=x[11]
K5=x[12]
V6=x[13]
K6=x[14]
k7=x[15]
K7=x[16]
k8=x[17]
K8=x[18]
V9=x[19]
K9=x[20]
V10=x[21]
K10=x[22]


sum((dMKKK-(V2*MKKK_P/(K2+MKKK_P)-V1*MKKK/((1+(MAPK_PP/Ki)^n)*(K1+MKKK))))^2+(dMKKK_P-(V1*MKKK/((1+(MAPK_PP/Ki)^n)*(K1+MKKK)) - V2*MKKK_P/(K2+MKKK_P)))^2+(dMKK-(V6*MKK_P/(K6+MKK_P) - k3*MKKK_P*MKK/(K3+MKK)))^2+(dMKK_P-(k3*MKKK_P*MKK/(K3+MKK) + V5*MKK_PP/(K5+MKK_PP) - k4*MKKK_P*MKK_P/(K4+MKK_P)  - V6*MKK_P/(K6+MKK_P)))^2+(dMKK_PP-(k4*MKKK_P*MKK_P/(K4+MKK_P)  - V5*MKK_PP/(K5+MKK_PP)))^2+(dMAPK-(V10*MAPK_P/(K10+MAPK_P)  - k7*MKK_PP*MAPK/(K7+MAPK)))^2+(dMAPK_P-(k7*MKK_PP*MAPK/(K7+MAPK) + V9*MAPK_PP/(K9+MAPK_PP) - k8*MKK_PP*MAPK_P/(K8+MAPK_P) - V10*MAPK_P/(K10+MAPK_P)))^2+(dMAPK_PP-(k8*MKK_PP*MAPK_P/(K8+MAPK_P)  - V9*MAPK_PP/(K9+MAPK_PP)))^2)

}


ui.m<-diag(1,22)
ci.v<-rep(1e-7,22)
P<-constrOptim(rep(5,22),min.function,NULL,ui=ui.m,ci=ci.v)

parameters<-P$par

names(parameters)<-c('V1','n','Ki','K1','V2','K2','k3','K3','k4','K4','V5','K5','V6','K6','k7','K7','k8','K8','V9','K9','V10','K10')

write.csv(data.matrix(parameters), file="set parameter solutions_egf path all points.csv")
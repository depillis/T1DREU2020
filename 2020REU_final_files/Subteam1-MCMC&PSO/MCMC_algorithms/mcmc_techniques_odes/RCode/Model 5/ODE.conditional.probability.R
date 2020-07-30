ode.conditional.probability<-function(ks,k.names,initial.data,observable.data,selected.species,forcing.list,dllname)
#ode.conditional.probability<-function(ks,required.data,alpha.var,beta.var,dllname)
{
library(deSolve)
parms<-ks[k.names]
#print(parms)
xstart<-initial.data[selected.species,]  #cambiado del 22 de julio de  2019
#print(xstart)
times<-sort(unique(observable.data[,'time']))
#print(xstart)

############## original solver######################3
#out1<-ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk_feedback", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsodes") ##ozi

###### function for problems with ode solver
##ODE23
##out1<-tryCatch(ode(xstart[selected.species],times,func = "derivs",parms = parms, dllname = dllname, initfunc = "initmod", nout = 1,outnames = "Sum",method="ode23"),error=function(e) -Inf)
##LSODEs
#out1<-tryCatch(ode(xstart[selected.species],times,func = "derivs",parms = parms, dllname = dllname, initfunc = "initmod", nout = 1,outnames = "Sum",method="lsodes"),error=function(e) -Inf)
 
   	##out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk_2005", initfunc = "initmod", nout = 1, outnames = "Sum",method="adams"),error=function(e) -Inf)
	###out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk_2005", initfunc = "initmod", nout = 1, outnames = "Sum",method="ode45"),error=function(e) -Inf)
	##out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk_2005", initfunc = "initmod", nout = 1, outnames = "Sum",method="ode23"),error=function(e) -Inf)
	##out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = "mapkk_2005", initfunc = "initmod", nout = 1, outnames = "Sum",method="lsoda"),error=function(e) -Inf)

#rk4,ode23,"ode45"
  Rec_i = ks['Rec_i'] 			#species
  xstart <- c(xstart[1],Rec_i,xstart[-1])
#scalar parameters
  scale_CD274mRNA_obs = ks['scale_CD274mRNA_obs']		#scaling parameter
  scale_IL13_cell_obs = ks['scale_IL13_cell_obs']		#scaling parameter
  scale_SOCS3mRNA_obs = ks['scale_SOCS3mRNA_obs']		#scaling parameter
  scale_pIL4Ra_obs = ks['scale_pIL4Ra_obs']				#scaling parameter
  scale_pJAK2_obs = ks['scale_pJAK2_obs']				#scaling parameter

#sigmaY1TimR <- ks['sigmaY1TimR',]  # it was computed from the experimental data, so it should not be sampled
#sigmaY2Step <- ks['sigmaY2Step',]  # it was computed from the experimental data, so it should not be sampled
#sigmaY2TimR <- ks['sigmaY2TimR',]  # it was computed from the experimental data, so it should not be sampled
#sigmaYDosR <- ks['sigmaYDosR',]  # it was computed from the experimental data, so it should not be sampled

simulated.data<-c()
exponential<-0

for(z in 1:4){
	forcing<-forcing.list[[z]]
	time.obs<-data.matrix( observable.data[observable.data[,1]==z,'time'])
	#    ode(xstart[selected.species],times,func = "derivs",parms = parms, dllname = dllname, initforc="forcc",forcings = forcing, initfunc = "initmod", nout = 1,outnames = "Sum")
	#	out1<-tryCatch(ode(xstart[selected.species],times,func = "derivs",parms = parms, dllname = dllname, initforc="forcc",forcings = forcing, initfunc = "initmod", nout = 1,outnames = "Sum"),error=function(e) -Inf)
	out1<-tryCatch(ode(xstart,times,func = "derivs",parms = parms, dllname = dllname, initforc="forcc",forcings = forcing, initfunc = "initmod", nout = 1,outnames = "Sum"),error=function(e) -Inf)
	
	#print(str(out1))
	if(length(out1)!=1){
		if(!all(times %in% out1[,1])){	exponential<--1e+308
		}else{
			temp.obs<-as.matrix(out1[out1[,'time']%in%time.obs,c('time','Rec','IL13_Rec','p_IL13_Rec','pJAK2','pSTAT5','SOCS3mRNA','SOCS3','CD274mRNA')])
			if(ncol(temp.obs)==1)temp.obs<-t(temp.obs)

			IL13_cell <- out1[out1[,'time']%in%time.obs,'IL13_Rec'] + out1[out1[,'time']%in%time.obs,'IL13_DecoyR'] + out1[out1[,'time']%in%time.obs,'p_IL13_Rec'] + out1[out1[,'time']%in%time.obs,'p_IL13_Rec_i']
			pIL4Ra <- out1[out1[,'time']%in%time.obs ,'p_IL13_Rec'] + out1[out1[,'time']%in%time.obs,'p_IL13_Rec_i']
                                                                                                                                                                                                                                           
			temp.obs<-cbind(temp.obs,IL13_cell, pIL4Ra)

			#colnames(observable.data)[c(2,seq(3,ncol(observable.data),3))]

			RecSurf_obs <-  temp.obs[,'Rec'] + temp.obs[,'IL13_Rec'] + temp.obs[,'p_IL13_Rec']
			IL13_cell_obs <-  scale_IL13_cell_obs * temp.obs[,'IL13_cell'] 
			pIL4Ra_obs <-  scale_pIL4Ra_obs * temp.obs[,'pIL4Ra'] 
			pJAK2_obs <-  scale_pJAK2_obs * temp.obs[,'pJAK2']
			SOCS3mRNA_obs <-  scale_SOCS3mRNA_obs * temp.obs[,'SOCS3mRNA'] 
			CD274mRNA_obs <-  scale_CD274mRNA_obs * temp.obs[,'CD274mRNA'] 
			SOCS3_obs <-  temp.obs[,'SOCS3'] 
			pSTAT5_obs <-  temp.obs[,'pSTAT5'] 

			temp.obs2<-cbind(z,temp.obs[,'time'],RecSurf_obs,IL13_cell_obs,pIL4Ra_obs,pJAK2_obs,SOCS3mRNA_obs,CD274mRNA_obs,SOCS3_obs,pSTAT5_obs)
			if(ncol(temp.obs2)==1){names(temp.obs2)[1:2]<-c('Model','time')
			}else{colnames(temp.obs2)[1:2]<-c('Model','time')}

			simulated.data<-rbind(simulated.data,temp.obs2)
		}
		}else{exponential<--1e+308}
}



rownames(simulated.data)<-NULL
if(exponential==0) {
	#IR1_P	IR1_P_std	IR1_P_ExpError	IRS1_P	IRS1_P_std	IRS1_P_ExpError
	obs.names <- colnames(observable.data)[seq(3,ncol(observable.data),3)]
	sd.names <- colnames(observable.data)[seq(4,ncol(observable.data),3)]
	fac.names <- colnames(observable.data)[seq(5,ncol(observable.data),3)]
	error.matrix<-sum(((observable.data[,obs.names]-simulated.data[,obs.names]*observable.data[,fac.names])^2)/(observable.data[,sd.names])^2)
	exponential<- -1*sum(error.matrix)
}
########################
##### Observables ######
#### Model 1 ####
# IR1_P<-k_IRP_1Step*(IRiP + IRp)
# IRS1_P<-IRSiP*k_IRSiP_1Step
#### Model 2 ####
# IRS1_P<-IRSiP*k_IRSiP_2Step
#### Model 3 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 4 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 5 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 6 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 7 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 8 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#### Model 9 ####
# IRS1_P<-IRSiP*k_IRSiP_DosR
#################
# parámetros de escalamiento
#define k_IRP_1Step parms[14]
#define k_IRSiP_1Step parms[15]
#define k_IRSiP_2Step parms[16]
#define k_IRSiP_DosR parms[17]
#define sigmaY1TimR parms[18]
#define sigmaY2Step parms[19]
#define sigmaY2TimR parms[20]
#define sigmaYDosR parms[21]

return(exponential)
}

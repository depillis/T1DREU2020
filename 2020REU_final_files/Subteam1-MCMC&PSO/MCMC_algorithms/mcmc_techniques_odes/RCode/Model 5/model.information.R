#source('model.information.R')
variance.vector<-(rep(5.6547672,23))^2  	#mass action
average.vector<-rep(-1.3697661,23)	#mass action
###  New hyperparameter computation for Michaelis Mente kinetics, mean=-2.62 & sd=3.46

magnitude.variance<-sqrt(sum(variance.vector*variance.vector))
#print(str(average.vector))



##########reading the kinetic parameters##########
Ks<-data.matrix(read.csv(paste(path.array[k,1],n.l,'/set parameter solutions_egf path steady state daniel.csv',sep=''),sep=',',row.names=1))
Ks<-Ks[c("CD274mRNA_production", "DecoyR_binding", "JAK2_p_inhibition", "JAK2_phosphorylation", "Kon_IL13Rec", "Rec_intern", "Rec_phosphorylation", "Rec_recycle", "SOCS3_accumulation", "SOCS3_degradation", "SOCS3_translation", "SOCS3mRNA_production", "STAT5_phosphorylation", "pJAK2_dephosphorylation", "pRec_degradation", "pRec_intern", "pSTAT5_dephosphorylation", "Rec_i", "scale_CD274mRNA_obs", "scale_IL13_cell_obs", "scale_SOCS3mRNA_obs", "scale_pIL4Ra_obs", "scale_pJAK2_obs"),1]
#Ks s a vector
fixed.par<-Ks
k.names<-c("CD274mRNA_production", "DecoyR_binding", "JAK2_p_inhibition", "JAK2_phosphorylation", "Kon_IL13Rec", "Rec_intern", "Rec_phosphorylation", "Rec_recycle", "SOCS3_accumulation", "SOCS3_degradation", "SOCS3_translation", "SOCS3mRNA_production", "STAT5_phosphorylation", "pJAK2_dephosphorylation", "pRec_degradation", "pRec_intern", "pSTAT5_dephosphorylation")
Ks.matrix<-t(array(Ks,dim=c(length(Ks),n.chains)))  #for metropolis
colnames(Ks.matrix) <- names(Ks)   ### for metropolis
log.current.sigma<-diag(0.1,length(Ks),length(Ks))  ### for adaptive
logik<-log(Ks)
###########################################################

################# reading experimental data #################
initial.values<-data.matrix(read.csv(paste(path.array[k,1],'/initial.species.values.csv',sep=''),sep=',',row.names=1))
observable.data<-data.matrix(read.csv(paste(path.array[k,1],'/Observable.data.csv',sep=''),sep=';'))
#selected.species<-c("Rec", "Rec_i", "IL13_Rec", "p_IL13_Rec", "p_IL13_Rec_i", "JAK2", "pJAK2", "STAT5", "pSTAT5", "SOCS3mRNA", "DecoyR", "IL13_DecoyR", "SOCS3", "CD274mRNA")
selected.species<-c("Rec", "IL13_Rec", "p_IL13_Rec", "p_IL13_Rec_i", "JAK2", "pJAK2", "STAT5", "pSTAT5", "SOCS3mRNA", "DecoyR", "IL13_DecoyR", "SOCS3", "CD274mRNA")
required.list<-list(initial.values,observable.data,selected.species)
#############################################################


######################################
## introduciendo los forcing values ##
######################################
times<-sort(unique(observable.data[,'time']))
forcing.matrix<-matrix(0,nrow=length(times),ncol=2)
forcing.matrix[,1]<-times
forcing.list<-list()
temp.list<-list()
temp.list[[1]]<-forcing.matrix
temp.list[[1]][,2] <- 91 ##SHP1	<- 91
forc.vector<-2.265*c(0, 20, 4, 80)  ## 2.265 is a constant value in the ODES that multiply the input il13_level
for(z in 1:4){
	temp.list[[2]]<-forcing.matrix
	temp.list[[2]][,2]<-forc.vector[z]
	forcing.list[[z]]<-temp.list
}


path.root<- "C:/"
my.folder<-'brenda_download/'  ## source --> https://www.brenda-enzymes.org/download_brenda_without_registration.php 
setwd(paste(path.root,my.folder,sep=''))

text.file<-readLines('brenda_download.txt')

text.matrix <- strsplit(text.file,'\t')
text.matrix.length <- sapply(1:length(text.matrix),function(x)length(text.matrix[[x]]))

depured.data<-text.file[text.matrix.length==2]

temp.data<-sapply(1:length(depured.data),function(x)unlist(strsplit(depured.data[x],'\t')))
temp.data<-t(temp.data)

########## Data to be conserved ###
#  KKM	Kcat/KM-Value substrate in {...}
#  KI	Ki-value	inhibitor in {...}
#  KM	KM-value	substrate in {...}
#  PR	protein
#################################
## I only selected information about species and parameters. Information regarding PH, EC number & others are excluded.

temp.data2<-temp.data[temp.data[,1]=="KKM"|temp.data[,1]=="KI"|
			temp.data[,1]=="KM"|temp.data[,1]=="PR",]


#[372835,] "KI" "#3# 12 {ADP}  (#3# 30ºC <13>) <13>"
#       "#1# Gallus gallus   <44>"
#[51,]  "#145# Lactococcus lactis subsp. lactis Q9CEN0 UniProt <250>"
#[180,]     "KM"   "#10# 50 {ethanol}  (#10# in presence of 50 mM SDS <192>) <192>"
#   <
# (
#  (

### matrix.split function for converting line data to matrix data
matrix.split<-function(z){
	if(temp.data2[z,1]=='PR'){  #z<-1
		temp.text<-unlist(strsplit(temp.data2[z,2],' (',fixed=TRUE))
		if(length(temp.text)==2) temp.text[2]<-paste('(',temp.text[2],sep='')
		if(length(temp.text)==1) temp.text<-unlist(strsplit(temp.data2[z,2],'   '))
		temp.text1<- c(temp.data2[z,1],unlist(strsplit(temp.text[1],'# ')),temp.text[2])  
	}
	if(temp.data2[z,1]=='KM'|temp.data2[z,1]=='KI'){ #z<-180
		temp.text<-unlist(strsplit(temp.data2[z,2]," {",fixed=TRUE))
		temp.text1<- c(temp.data2[z,1], unlist(strsplit(temp.text[1],'# ')),paste('{',temp.text[2],sep='')) 
	}
	if(temp.data2[z,1]=='KKM'){ #z<-180
		temp.text<-unlist(strsplit(temp.data2[z,2]," (",fixed=TRUE))
		temp.text1<- c(temp.data2[z,1], unlist(strsplit(temp.text[1],'# ')),paste('(',temp.text[2],sep='')) 
		temp.text2<-unlist(strsplit(temp.text1[3],' <'))
		if(length(temp.text2)==2) temp.text1[3:4] <- c(temp.text1[1],paste('<',temp.text2[2],' ',temp.text1[4],sep=''))
	}
	return(unlist(temp.text1))
}


temp.data3<-sapply(1:nrow(temp.data2), matrix.split )
previous.temp.data2<-temp.data2
amount.data <- sapply(1:length(temp.data3), function(x) length(temp.data3[[x]]))
temp.data2<-temp.data2[amount.data==4,]
temp.data4<-sapply(1:nrow(temp.data2), matrix.split )
temp.data4<-t(temp.data4)

## for knowinng where the parameters line starts and end ##
k.rows<-which(temp.data4[,1]!='PR')
rows.difference<-k.rows[-1]-k.rows[-length(k.rows)]
k.rows.1<-k.rows[c(1,which(rows.difference!=1)+1)]
## for knowinng where the proteins' line startas and end ###
pr.rows<-which(temp.data4[,1]=='PR')
pr.rows.difference<-pr.rows[-1]-pr.rows[-length(pr.rows)]
pr.rows.1<-pr.rows[c(1,which(pr.rows.difference!=1)+1)]


### loops for substracting only parameters cells and adding a extra column with the protein name
new.matrix<-c()
for(i in 1:(length(k.rows.1)-1)){
	rows.vector   <-k.rows.1[i]:(pr.rows.1[i+1]-1)
	pr.rows.vector<-pr.rows.1[i]:(  k.rows.1[i]-1)
	ID.matrix <- temp.data4[pr.rows.vector,]
	if(length(ID.matrix)==4) ID.matrix <- t(as.matrix(ID.matrix))
	for(j in rows.vector){
		ID.vector<- unlist(strsplit(temp.data4[j,2],','))
		ID<-ID.vector[1]
		new.matrix<-rbind(new.matrix,c(temp.data4[j,],ID.matrix[ID.matrix[,2]==ID,3]))
		if(length(ID.vector)!=1) for(ID in ID.vector[-1])new.matrix<-rbind(new.matrix,c(temp.data4[j,],ID.matrix[ID.matrix[,2]==paste('#',ID,sep=''),3]))
	}
}


new.matrix3<-new.matrix[new.matrix[,3]!='KKM',]  #cleanning lines where there are shown the characters 'KKM' instead of k values

###for indetifying a k range values intesd of k values
range.index<-sapply(1:length(new.matrix3[,3]),function(x) is.na(as.numeric(new.matrix3[x,3])) )
range.vector<-which(range.index==TRUE)

############## converting the range values in lines with just values
temp.matrix<-new.matrix3
for(i.range in range.vector)new.matrix3<-rbind(new.matrix3,cbind( rbind(new.matrix3[i.range,1:2],new.matrix3[i.range,1:2]),unlist(strsplit(new.matrix3[i.range,3],'-')), rbind(new.matrix3[i.range,4:5], new.matrix3[i.range,4:5])))
new.matrix3<-new.matrix3[-range.vector,]
###########

################# Identifying rows with homo sapiens data #############3
homo.sapiens.vector<- sapply(1:length(new.matrix3[,3]),function(x)  grepl('omo sapien',new.matrix3[x,5]))  #for recognizing "omo sapien" pattern
homo.sapiens.matrix<-new.matrix3[homo.sapiens.vector,]
################################################################3



library(MASS)
KKM.mean<-fitdistr(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KKM',3])),'normal')[[1]][1]
KKM.sd<-fitdistr(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KKM',3])),'normal')[[1]][2]
KM.mean<-fitdistr(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KM',3])),'normal')[[1]][1]
KM.sd<-fitdistr(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KM',3])),'normal')[[1]][2]
KI.mean<-fitdistr(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KI',3])),'normal')[[1]][1]
KI.sd<-fitdistr(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KI',3])),'normal')[[1]][2]

prior.data<-rbind(c(KKM.mean,KKM.sd),c(KM.mean,KM.sd),c(KI.mean,KI.sd))
colnames(prior.data)<-c('MEAN','SD')
rownames(prior.data)<-c('KKM','KM','KI')

write.csv(prior.data,file='Brenda hyperparameters values for kinetics parameter_Homo sapiens.csv')

pdf('histogram from brenda.pdf')
#d<-density(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KKM',3])))
x<-unique(sort(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KKM',3]))))
y<-(exp(-(unique(x)-KKM.mean)^2/(2*pi^2)))/sqrt(2*pi*KKM.sd^2) #density function
hist(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KKM',3])),main=paste('Ln KKM ~ (',round(KKM.mean,2),',',round(KKM.sd,2),')',sep=''),ylab='Frequency',xlab='Ln KMM values',col='lightblue')
hist(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KKM',3])),main=paste('Ln KKM ~ (',round(KKM.mean,2),',',round(KKM.sd,2),')',sep=''),ylab='Frequency',xlab='Ln KMM values',col='pink',freq=FALSE)
lines(x,y)

x<-unique(sort(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KM',3]))))
y<-(exp(-(unique(x)-KM.mean)^2/(2*pi^2)))/sqrt(2*pi*KM.sd^2) #density function
hist(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KM',3])),main=paste('Ln KM ~ (',round(KM.mean,2),',',round(KM.sd,2),')',sep=''),ylab='Frequency',xlab='Ln KM values',col='lightblue')
hist(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KM',3])),main=paste('Ln KM ~ (',round(KM.mean,2),',',round(KM.sd,2),')',sep=''),ylab='Frequency',xlab='Ln KM values',col='pink',freq=FALSE)
lines(x,y)

x<-unique(sort(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KI',3]))))
y<-(exp(-(unique(x)-KI.mean)^2/(2*pi^2)))/sqrt(2*pi*KI.sd^2) #density function
hist(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KI',3])),main=paste('Ln KI ~ (',round(KI.mean,2),',',round(KI.sd,2),')',sep=''),ylab='Frequency',xlab='Ln KI values',col='lightblue')
hist(log(as.numeric(homo.sapiens.matrix[as.numeric(homo.sapiens.matrix[,3])>=0 & homo.sapiens.matrix[,1]=='KI',3])),main=paste('Ln KI ~ (',round(KI.mean,2),',',round(KI.sd,2),')',sep=''),ylab='Frequency',xlab='Ln KI values',col='pink',freq=FALSE)
lines(x,y)
dev.off()



#source('data.format.R')

proteins<-data.matrix(read.csv("noise data 0.1%.csv", sep="\t"))
#proteins<-data.matrix(read.csv("noise data 0.1%.csv", row.names=1))
initial.proteins<-proteins[proteins[,'time']==0,]

proteins1<-data.matrix(read.csv("noise data 0.1%_r1.csv", row.names=1))
proteins2<-data.matrix(read.csv("noise data 0.1%_r2.csv", row.names=1))
proteins3<-data.matrix(read.csv("noise data 0.1%_r3.csv", row.names=1))

sigma.data<-data.matrix(read.csv("variance noise data 0.1%.csv", row.names=1))

prior.solution.steady<-data.matrix(read.csv("set parameter solutions_egf path steady state daniel.csv", row.names=1,sep="\t"))
prior.solution.all<-data.matrix(read.csv("set parameter solutions_egf path all points.csv", row.names=1,sep="\t"))

all.data<-list(initial.proteins,proteins1,proteins2,proteins3,sigma.data,prior.solution.steady,prior.solution.all)
### [[1]] vector, [[2]] matrix n.timexn.species, [[3]] matrix n.timexn.species, [[4]] matrix n.timexn.species, [[5]] matrix n.timexn.species, [[6]] matrix n.parameterx1,[[7]] matrix n.parameterx1
save(all.data,file='input.data.a.RData')






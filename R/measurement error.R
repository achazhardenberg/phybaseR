### Variability in trait model

# As an example of a model taking account of data for which repeated measures
# are available for each species we add white noise to the rhino data.
# To do this each trait will be transformed to a distribution with mean equal
# to the original species trait value and SD equal to a proportion (25%) of the among-species variance for that trait.

##Simulation of repeated measures for each trait
set.seed(12345)

rhino.dat<-read.csv("rhino.csv")

BM.new<-data.frame()
for(i in 1:length(rhino.dat$BM)){
  BM.mult<-sample(rnorm(500,mean=rhino.dat$BM[i],sd=sd(rhino.dat$BM)/4),10)
  BM.new<-rbind(BM.new,BM.mult)
}
names(BM.new)<-c("BM1","BM2","BM3","BM4","BM5","BM6","BM7","BM8","BM9","BM10")

NL.new<-data.frame()
for(i in 1:length(rhino.dat$NL)){
  NL.mult<-sample(rnorm(500,mean=rhino.dat$NL[i],sd=sd(rhino.dat$NL)/4),10)
  NL.new<-rbind(NL.new,NL.mult)
}
names(NL.new)<-c("NL1","NL2","NL3","NL4","NL5","NL6","NL7","NL8","NL9","NL10")

LS.new<-data.frame()
for(i in 1:length(rhino.dat$LS)){
  LS.mult<-sample(rnorm(500,mean=rhino.dat$LS[i],sd=sd(rhino.dat$LS)/4),10)
  LS.new<-rbind(LS.new,LS.mult)
}
names(LS.new)<-c("LS1","LS2","LS3","LS4","LS5","LS6","LS7","LS8","LS9","LS10")

DD.new<-data.frame()
for(i in 1:length(rhino.dat$DD)){
  DD.mult<-sample(rnorm(500,mean=rhino.dat$DD[i],sd=sd(rhino.dat$DD)/4),10)
  DD.new<-rbind(DD.new,DD.mult)
}
names(DD.new)<-c("DD1","DD2","DD3","DD4","DD5","DD6","DD7","DD8","DD9","DD10")

RS.new<-data.frame()
for(i in 1:length(rhino.dat$RS)){
  RS.mult<-sample(rnorm(500,mean=rhino.dat$RS[i],sd=sd(rhino.dat$RS)/4),10)
  RS.new<-rbind(RS.new,RS.mult)
}
names(RS.new)<-c("RS1","RS2","RS3","RS4","RS5","RS6","RS7","RS8","RS9","RS10")

data<-cbind(BM.new,NL.new,LS.new,DD.new,RS.new)

BM<-subset(data, select=BM1:BM10)
BM<-data.matrix(BM, rownames.force=NA)
BM2<-t(BM)
BMMulti<-as.vector(c(BM2[,1:100]))

NL<-subset(data, select=NL1:NL10)
NL<-data.matrix(NL, rownames.force=NA)
NL2<-t(NL)
NLMulti<-as.vector(c(NL2[,1:100]))

LS<-subset(data, select=LS1:LS10)
LS<-data.matrix(LS, rownames.force=NA)
LS2<-t(LS)
LSMulti<-as.vector(c(LS2[,1:100]))

DD<-subset(data, select=DD1:DD10)
DD<-data.matrix(DD, rownames.force=NA)
DD2<-t(DD)
DDMulti<-as.vector(c(DD2[,1:100]))

RS<-subset(data, select=RS1:RS10)
RS<-data.matrix(RS, rownames.force=NA)
RS2<-t(RS)
RSMulti<-as.vector(c(RS2[,1:100]))

dataMulti<-as.data.frame(cbind(BMMulti,NLMulti,LSMulti,DDMulti,RSMulti))

SP<-paste("s",sort(rep(1:100,10)),sep="")

dataMulti<-cbind(SP,dataMulti)
RhinoMulti.dat <- dataMulti


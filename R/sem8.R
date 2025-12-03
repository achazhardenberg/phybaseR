library(ape)
library(caper)
library(nlme)
library(rjags)
library(R2jags)

rhino.dat<-read.csv("rhino.csv")
rhino.tree<-read.tree("rhino.tree")

rhino.tree$edge.length<-rhino.tree$edge.length/max(branching.times(rhino.tree))

rhino.vcv <- vcv.phylo(rhino.tree)
ID<-diag(100)
VCV<-rhino.vcv

nc<-3  #number of chains
ni<-12000 #number of iterations
nb<-2000 #number of burnin
nt<-10 #thinning number

set.seed(12345) #we set a seed so results from the mcmc are perfectly replicable

sem.data<-list(BM=rhino.dat$BM,NL=rhino.dat$NL,DD=rhino.dat$DD,LS=rhino.dat$LS, RS=rhino.dat$RS,VCV=VCV,ID=ID,Nspec=100)

sem8.jg<-function(){
  #Structural equations 
  for (i in 1:Nspec) {
    muLS[i] <- alphaLS+betaBM*BM[i]
    muNL[i] <- alphaNL+betaBM2*BM[i]+betaRS*RS[i]
    muDD[i] <- alphaDD+betaNL*NL[i]
  }
  #Multivariate normal likelihoods
  LS[1:Nspec] ~ dmnorm(muLS[],TAUls)
  NL[1:Nspec] ~ dmnorm(muNL[],TAUnl)
  DD[1:Nspec] ~ dmnorm(muDD[],TAUdd)
  #Priors
  alphaLS ~ dnorm(0,1.0E-06)
  alphaNL ~ dnorm(0,1.0E-06)
  alphaDD ~ dnorm(0,1.0E-06)
  betaBM ~ dnorm(0,1.0E-06)
  betaBM2 ~ dnorm(0,1.0E-06)
  betaRS ~ dnorm(0,1.0E-06)
  betaNL ~ dnorm(0,1.0E-06)
  lambdaLS ~ dunif(0,1)
  lambdaNL ~ dunif(0,1)
  lambdaDD ~ dunif(0,1)
  tauLS ~ dgamma(1,1)
  tauNL ~ dgamma(1,1)
  tauDD ~ dgamma(1,1)
  sigmaLS <- 1/sqrt(tauLS)
  sigmaNL <- 1/sqrt(tauNL)
  sigmaDD <- 1/sqrt(tauDD)
  #lambda computation
  MlamLS <- lambdaLS*VCV+(1-lambdaLS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  MlamNL <- lambdaNL*VCV+(1-lambdaNL)*ID
  TAUnl <- tauNL*inverse(MlamNL)
  MlamDD <- lambdaDD*VCV+(1-lambdaDD)*ID
  TAUdd <- tauDD*inverse(MlamDD)
}

params8 <- c("alphaLS", "alphaNL","alphaDD","betaBM","betaBM2","betaRS","betaNL")
sem8.mcmc<-jags.parallel(data=sem.data,model.file=sem8.jg,n.chains=3,n.iter=12000,n.burnin=2000, n.thin=10,parameters.to.save=params8)

samples.sem8 <- jags.samples(sem8.mcmc$model, 
                             c("WAIC","deviance"), 
                             type = "mean", 
                             n.iter = ni,
                             n.burnin = nb,
                             n.thin = nt)

samples.sem8$p_waic <- samples.sem8$WAIC
samples.sem8$waic <- samples.sem8$deviance + samples.sem8$p_waic
tmp <- sapply(samples.sem8, sum)
waic.sem8 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.sem8

save(file="sem8.Rdata", list=c("sem8.mcmc","samples.sem8","waic.sem8"))

# This file shows how to get simulated data set named "simulation2.RData".

# Now we show how to simulate data for the second application in the paper Application 2: Insurance Company Loss and Expense (Reinsurance)

# There are two copulas applied. Assuming a combined deductible is applied to the sum of loss and ALAE.
library(copula)
library(VGAM)
set.seed(1211)
nsim<-5000000 # The number of simulation. Feel free to change it. But keep in mind that a lower number of simulation may give you inaccurate result in the end.
X_G <- matrix(0,nsim,2)-> X_iid-> X_F->X_comon->X_counter
SumXYlA <- matrix(0,nsim,3) ->SumXYlA_f

#independent coppula
copulaParamid<-iRho(normalCopula(param=1, dim = 2), rho=0)
ncid<-normalCopula(param=copulaParamid, dim =2)
#nccomon<-normalCopula(param=1, dim =2) # comonotonic copula
#nccounter<-normalCopula(param=-1, dim =2) # countermonotonic copula

#iid severity of loss and ALAE
Xid<-rCopula(nsim, ncid)
X_iid[,1] <- qparetoII(Xid[,1], scale = 14453, shape=1.135) #severity of loss
X_iid[,2] <- qparetoII(Xid[,2], scale = 15133, shape=2.223) #severity of ALAE

# comonotonic version severity of loss and ALAE, used in supplentary material
#Xcomon<-rCopula(nsim, nccomon)
#X_comon[,1] <- qparetoII(Xcomon[,1], scale = 14453, shape=1.135) #severity of loss
#X_comon[,2] <- qparetoII(Xcomon[,2], scale = 15133, shape=2.223) #severity of ALAE

# countermonotonic version severity of loss and ALAE, used in supplentary material
#Xcounter<-rCopula(nsim, nccounter)
#X_counter[,1] <- qparetoII(Xcounter[,1], scale = 14453, shape=1.135) #severity of loss
#X_counter[,2] <- qparetoII(Xcounter[,2], scale = 15133, shape=2.223) #severity of ALAE

#Gumbel copula
gumbel.cop <- gumbelCopula(param=1.453,dim=2) # the parameters of Gumbel copula is from Jed's paper; see "Understanding relationships using copulas".
rho(gumbel.cop)
GumbelParam <- iRho(gumbelCopula(dim = 2), rho=rho(gumbel.cop))
ncg <- gumbelCopula(param = GumbelParam, dim = 2) 
Xg <-rCopula(nsim, ncg)
X_G[,1] <- qparetoII(Xg[,1], scale = 14036, shape=1.112) #severity of loss from Gumbel copula
X_G[,2] <- qparetoII(Xg[,2], scale = 14219, shape=2.118) #severity of ALAE from Gumbel copula

#Frank copula
frank.cop <- frankCopula(param=3.158,dim=2) # the parameters of Frank copula is from Jed's paper; see "Understanding relationships using copulas".
FrankParam <- iRho(frankCopula(param=1, dim = 2), rho=rho(frank.cop))
ncf <- frankCopula(param = FrankParam, dim = 2)
Xf <-rCopula(nsim, ncf)
X_F[,1] <- qparetoII(Xf[,1], scale = 14558, shape=1.115) #severity of loss
X_F[,2] <- qparetoII(Xf[,2], scale = 16678, shape=2.309) #severity of ALAE

SumXYlA[,1]<-X_iid[,1]+X_iid[,2] #iid sum of two severities
SumXYlA[,2]<-X_G[,1]+X_G[,2] #gumbel sum of two severities
SumXYlA[,3]<-X_F[,1]+X_F[,2] #frank sum of two severities
#Sumcomon<-X_comon[,1]+X_comon[,2] #comonotonic version sum of two severities
#Sumcounter<-X_counter[,1]+X_counter[,2] #countermonotonic version sum of two severities


dlA<- 0:100*500 # combined deductibles 0 to 50000
dlA[1]<-0.000001
dlA[2]<-0.0001
dlA[3]<-0.1
LlA<-matrix(0,length(dlA),3)->LlAp
Lcomon<-matrix(0,length(dlA),1)->Lcounter

for (i in 1:3){
  for (j in 1:length(dlA)) {
    d = dlA[j]
    LlA[j,i]  <- mean(pmin(SumXYlA[,i],d)) # cost of deductible for iid, gumbel, and frank cases
    LlAp[j,i]  <-mean(SumXYlA[,i])-mean(pmin(SumXYlA[,i],d)) # premium  for iid, gumbel, and frank cases
  }
}

#comon & counter version
#for (j in 1:length(dlA)) {
#  d = dlA[j]
#  Lcomon[j,1]  <- mean(pmin(Sumcomon,d)) # cost of deductible comon
#  Lcounter[j,1]  <-mean(pmin(Sumcounter,d))  # cost of deductible counter
#}


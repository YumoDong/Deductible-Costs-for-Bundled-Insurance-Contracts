# This file shows how to get simulated data set named "simulation3.RData".
# Now we show how to simulate data for the third application in the paper Application 3: Quotes of a Pet Insurance Contract (Personal)
library(optimx)
library(nlme) 
library(VGAM)
library(copula)
library(statmod)
library(numDeriv)
library(DEoptim)

# The first part of this file is to show how to use DEoptim() to estimate the parameters of assumed distributions accoding to MAPE; see Table 7 of the paper

################################## First case: gamma + gamma; see Table 7 Model 1.############################
set.seed(2020202)
nite<-200 #number of Interation is 200 for DEoptim
f = function(par){
  # quotes data of pet insurance with combined limits, combined deductibles, and reimbursements.
  cl<-c(5000,8000,10000,15000) #4 combined limits
  cd<-c(200,300,500,750,1000)  #5 combined deductibles
  price_real0.9<-matrix(c(847.2,	974.4,	1059,	1270.92,
                          728.64,	837.84,	910.8,	1092.96,
                          550.68,	633.24,	688.44,	826.08,
                          355.8,	409.2,	444.84,	533.76,
                          254.16,	292.32,	317.76,	381.24),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.9
  
  price_real0.8<-matrix(c(627.6,	721.68,	784.56,	941.4,
                          539.76,	620.76,	674.64,	809.64,
                          407.88,	469.08,	509.88,	612,
                          263.52,	303.12,	329.4, 395.4,
                          188.28,	216.48,	235.32,	282.36),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.8
  
  price_real0.7<-matrix(c(546, 627.96,	682.56,	819,
                          469.56,	540,	587.04,	704.4,
                          354.84,	408.12,	443.64,	532.32,
                          229.32,	263.64,	286.68,	343.92,
                          163.8,	188.28,	204.72,	245.64),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.7
  
  remrate<-c(0.9, 0.8, 0.7) #reimbursement rate
  mse<-rep(0,4)
  # Assumptions
  nsim<-500000 #number of simulation
  set.seed(2020202)
  price_estimated0.9<- matrix(0,length(cd),length(cl))->price_estimated0.8->price_estimated0.7
  X_md   <- matrix(0,nsim,2) ->freq_md
  
  copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
  nc <- normalCopula(param=copulaParam, dim =2)
  X <- rCopula(nsim, nc)
  Xf<- rCopula(nsim, normalCopula(param=iRho(normalCopula(param=1, dim = 2), rho=0.8), dim =2)) #rho of frequency coopula is fixed at 0.8, so we don't estimate it.
  
  X_md[,1] <- qgamma(X[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1
  X_md[,2] <- qgamma(X[,2],scale = par[3]/par[1], shape=par[1]) #severity of X2
  freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1
  freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2
  SumXY <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim
  for (i in 1:length(cd)) {
    for (j in 1:length(cl)){
      d = cd[i]
      l = cl[j]
      loss0.9_XY <- pmin(pmax(SumXY-d,0)*remrate[1] ,l)
      loss0.8_XY <- pmin(pmax(SumXY-d,0)*remrate[2] ,l)
      loss0.7_XY <- pmin(pmax(SumXY-d,0)*remrate[3] ,l)
      price_estimated0.9[i,j]<-mean(loss0.9_XY)
      price_estimated0.8[i,j]<-mean(loss0.8_XY)
      price_estimated0.7[i,j]<-mean(loss0.7_XY)
    }
  }
  mape<- sum(sqrt(((price_estimated0.9-price_real0.9)/price_real0.9)^2) + sqrt(((price_estimated0.8-price_real0.8)/price_real0.8)^2) +sqrt(((price_estimated0.7-price_real0.7)/price_real0.7)^2))/(length(cd)*length(cl)*length(remrate))
  
  mape #mape
}

Mapfun <- mappingFun <- function(x) { # this function defines the decimal place of each parameter
  x[1] <- round(x[1], 2) #shape
  x[2] <- round(x[2], 2) #shape
  x[3] <- round(x[3], 1) #mean
  x[4] <- round(x[4], 1) #mean
  x[5] <- round(x[5], 2) #rho of copula
  x
}

result1=DEoptim(fn = f, upper = c(20,20,3000,3000,0.95), lower = c(0.5,0.5,500,500,0.4),fnMap = Mapfun, DEoptim.control(itermax = nite, packages = c("copula","VGAM"), parallelType=1)) #parallel mode is used to speed up
# after 200 interations, you can see esimated shape (accident), shape (illness), mean (accident), mean (illness), and rho of copula.

################################## Second case: gamma + Pareto; see Table 7 Model 2############################
f2 = function(par){
  cl<-c(5000,8000,10000,15000) #4 combined limits
  cd<-c(200,300,500,750,1000)  #5 combined deductible
  price_real0.9<-matrix(c(847.2,	974.4,	1059,	1270.92,
                          728.64,	837.84,	910.8,	1092.96,
                          550.68,	633.24,	688.44,	826.08,
                          355.8,	409.2,	444.84,	533.76,
                          254.16,	292.32,	317.76,	381.24),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.9
  
  price_real0.8<-matrix(c(627.6,	721.68,	784.56,	941.4,
                          539.76,	620.76,	674.64,	809.64,
                          407.88,	469.08,	509.88,	612,
                          263.52,	303.12,	329.4, 395.4,
                          188.28,	216.48,	235.32,	282.36),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.8
  
  price_real0.7<-matrix(c(546, 627.96,	682.56,	819,
                          469.56,	540,	587.04,	704.4,
                          354.84,	408.12,	443.64,	532.32,
                          229.32,	263.64,	286.68,	343.92,
                          163.8,	188.28,	204.72,	245.64),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.7
  
  remrate<-c(0.9, 0.8, 0.7) #reimbursement rate
  mse<-rep(0,4)
  # Assumptions
  nsim<-500000#number of simulation
  set.seed(2020202)
  price_estimated0.9<- matrix(0,length(cd),length(cl))->price_estimated0.8->price_estimated0.7
  X_md   <- matrix(0,nsim,2)->freq_md
  
  copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
  nc <- normalCopula(param=copulaParam, dim =2)
  X <- rCopula(nsim, nc)
  Xf<- rCopula(nsim, normalCopula(param=iRho(normalCopula(param=1, dim = 2), rho=0.8), dim =2))
  
  X_md[,1] <- qgamma(X[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1
  # Weibull(shape = a , scale = s); E(X) = s*(1/a)!
  X_md[,2] <- qparetoII(X[,2],scale = par[3]*(par[1]-1), shape=par[1]) #severity of X2
  freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1
  freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2
  SumXY <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim
  
  for (i in 1:length(cd)) {
    for (j in 1:length(cl)){
      d = cd[i]
      l = cl[j]
      loss0.9_XY <- pmin(pmax(SumXY-d,0)*remrate[1] ,l)
      loss0.8_XY <- pmin(pmax(SumXY-d,0)*remrate[2] ,l)
      loss0.7_XY <- pmin(pmax(SumXY-d,0)*remrate[3] ,l)
      price_estimated0.9[i,j]<-mean(loss0.9_XY)
      price_estimated0.8[i,j]<-mean(loss0.8_XY)
      price_estimated0.7[i,j]<-mean(loss0.7_XY)
    }
  }
  
  mape<- sum(sqrt(((price_estimated0.9-price_real0.9)/price_real0.9)^2) + sqrt(((price_estimated0.8-price_real0.8)/price_real0.8)^2) +sqrt(((price_estimated0.7-price_real0.7)/price_real0.7)^2))/(length(cd)*length(cl)*length(remrate))
  
  mape
}

result2=DEoptim(fn = f2, upper = c(20,20,3000,3000,1), lower = c(2.01,0.5,500,500,0.4),fnMap = Mapfun, DEoptim.control(itermax = nite, packages = c("copula","VGAM"), parallelType=1)) #parallel mode is used to speed up
# after 200 interations, you can see esimated shape (accident), shape (illness), mean (accident), mean (illness), and rho of copula.

################################## Third case: gamma + Weibull; see Table 7 Model 3############################
f3 = function(par){
  cl<-c(5000,8000,10000,15000) #4 combined limits
  cd<-c(200,300,500,750,1000)  #5 combined deductible
  price_real0.9<-matrix(c(847.2,	974.4,	1059,	1270.92,
                          728.64,	837.84,	910.8,	1092.96,
                          550.68,	633.24,	688.44,	826.08,
                          355.8,	409.2,	444.84,	533.76,
                          254.16,	292.32,	317.76,	381.24),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.9
  
  price_real0.8<-matrix(c(627.6,	721.68,	784.56,	941.4,
                          539.76,	620.76,	674.64,	809.64,
                          407.88,	469.08,	509.88,	612,
                          263.52,	303.12,	329.4, 395.4,
                          188.28,	216.48,	235.32,	282.36),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.8
  
  price_real0.7<-matrix(c(546, 627.96,	682.56,	819,
                          469.56,	540,	587.04,	704.4,
                          354.84,	408.12,	443.64,	532.32,
                          229.32,	263.64,	286.68,	343.92,
                          163.8,	188.28,	204.72,	245.64),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.7
  
  remrate<-c(0.9, 0.8, 0.7) #reimbursement rate
  mse<-rep(0,4)
  # Assumptions
  
  nsim<-500000#number of simulation
  set.seed(2020202)
  price_estimated0.9<- matrix(0,length(cd),length(cl))->price_estimated0.8->price_estimated0.7
  X_md   <- matrix(0,nsim,2)->freq_md
  copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
  nc <- normalCopula(param=copulaParam, dim =2)
  X <- rCopula(nsim, nc)
  Xf<- rCopula(nsim, normalCopula(param=iRho(normalCopula(param=1, dim = 2), rho=0.8), dim =2))
  
  X_md[,1] <- qgamma(X[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1
  # Weibull(shape = a , scale = s); E(X) = s*(1/a)!
  X_md[,2] <- qweibull(X[,2],scale = par[3]/gamma(1/par[1]+1), shape=par[1]) #severity of X2
  freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1
  freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2
  SumXY <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim
  
  
  for (i in 1:length(cd)) {
    for (j in 1:length(cl)){
      d = cd[i]
      l = cl[j]
      loss0.9_XY <- pmin(pmax(SumXY-d,0)*remrate[1] ,l)
      loss0.8_XY <- pmin(pmax(SumXY-d,0)*remrate[2] ,l)
      loss0.7_XY <- pmin(pmax(SumXY-d,0)*remrate[3] ,l)
      price_estimated0.9[i,j]<-mean(loss0.9_XY)
      price_estimated0.8[i,j]<-mean(loss0.8_XY)
      price_estimated0.7[i,j]<-mean(loss0.7_XY)
    }
  }
  mape<- sum(sqrt(((price_estimated0.9-price_real0.9)/price_real0.9)^2) + sqrt(((price_estimated0.8-price_real0.8)/price_real0.8)^2) +sqrt(((price_estimated0.7-price_real0.7)/price_real0.7)^2))/(length(cd)*length(cl)*length(remrate))
  
  mape
}

result3=DEoptim(fn = f3, upper = c(10,15,3000,3000,0.95), lower = c(0.1,0.5,500,500,0.4),fnMap = Mapfun, DEoptim.control(itermax = nite, packages = c("copula","VGAM"), parallelType=1)) #parallel mode is used to speed up
# after 200 interations, you can see esimated shape (accident), shape (illness), mean (accident), mean (illness), and rho of copula.

################################## Third case: Weibull + Weibull; see Table 7 Model 4############################
f4 = function(par){
  cl<-c(5000,8000,10000,15000) #4 combined limits
  cd<-c(200,300,500,750,1000)  #5 combined deductible
  price_real0.9<-matrix(c(847.2,	974.4,	1059,	1270.92,
                          728.64,	837.84,	910.8,	1092.96,
                          550.68,	633.24,	688.44,	826.08,
                          355.8,	409.2,	444.84,	533.76,
                          254.16,	292.32,	317.76,	381.24),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.9
  
  price_real0.8<-matrix(c(627.6,	721.68,	784.56,	941.4,
                          539.76,	620.76,	674.64,	809.64,
                          407.88,	469.08,	509.88,	612,
                          263.52,	303.12,	329.4, 395.4,
                          188.28,	216.48,	235.32,	282.36),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.8
  
  price_real0.7<-matrix(c(546, 627.96,	682.56,	819,
                          469.56,	540,	587.04,	704.4,
                          354.84,	408.12,	443.64,	532.32,
                          229.32,	263.64,	286.68,	343.92,
                          163.8,	188.28,	204.72,	245.64),nrow=5,ncol=4,byrow = TRUE) # Annual premium for Reimbursement=0.7
  
  remrate<-c(0.9, 0.8, 0.7) #reimbursement rate
  mse<-rep(0,4)
  # Assumptions
  nsim<-500000#number of simulation
  set.seed(2020202)
  price_estimated0.9<- matrix(0,length(cd),length(cl))->price_estimated0.8->price_estimated0.7
  X_md   <- matrix(0,nsim,2)->freq_md
  copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
  nc <- normalCopula(param=copulaParam, dim =2)
  X <- rCopula(nsim, nc)
  Xf<- rCopula(nsim, normalCopula(param=iRho(normalCopula(param=1, dim = 2), rho=0.8), dim =2))
  
  X_md[,1] <- qweibull(X[,1],scale = par[4]/gamma(1/par[2]+1), shape=par[2])  #severity of X1
  # Weibull(shape = a , scale = s); E(X) = s*(1/a)!
  X_md[,2] <- qweibull(X[,2],scale = par[3]/gamma(1/par[1]+1), shape=par[1]) #severity of X2
  freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1
  freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2
  SumXY <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim
  
  for (i in 1:length(cd)) {
    for (j in 1:length(cl)){
      d = cd[i]
      l = cl[j]
      loss0.9_XY <- pmin(pmax(SumXY-d,0)*remrate[1] ,l)
      loss0.8_XY <- pmin(pmax(SumXY-d,0)*remrate[2] ,l)
      loss0.7_XY <- pmin(pmax(SumXY-d,0)*remrate[3] ,l)
      price_estimated0.9[i,j]<-mean(loss0.9_XY)
      price_estimated0.8[i,j]<-mean(loss0.8_XY)
      price_estimated0.7[i,j]<-mean(loss0.7_XY)
    }
  }
  mape<- sum(sqrt(((price_estimated0.9-price_real0.9)/price_real0.9)^2) + sqrt(((price_estimated0.8-price_real0.8)/price_real0.8)^2) +sqrt(((price_estimated0.7-price_real0.7)/price_real0.7)^2))/(length(cd)*length(cl)*length(remrate))
  
  mape
}

result4=DEoptim(fn = f4, upper= c(15,15,3000,3000,0.95), lower = c(0.1,0.1,500,500,0.4),fnMap = Mapfun, DEoptim.control(itermax = nite, packages = c("copula","VGAM"), parallelType=1)) #parallel mode is used to speed up
# after 200 interations, you can see esimated shape (accident), shape (illness), mean (accident), mean (illness), and rho of copula.


# The second part of this file is to show how to use estimated the parameters to get simulated data set called "simulation3.RData".

set.seed(1211)
nsim<-2500000 # number of simulation is 2500000
#X1 illness
#X2 accident
X_md   <- matrix(0,nsim,2)->freq_md->X_mdid->freq_mdid->X_common->freq_common->X_counter->freq_counter
SumXY <- matrix(0,nsim,4)->SumXYid
#copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
#nc <- normalCopula(param=copulaParam, dim =2)
copulaParamid<-iRho(normalCopula(param=1, dim = 2), rho=0)
ncid<-normalCopula(param=copulaParamid, dim =2)
Xf<- rCopula(nsim, normalCopula(param=iRho(normalCopula(param=1, dim = 2), rho=0.8), dim =2)) #frequency normal copula 0.8 rho
Xid<-rCopula(nsim, ncid) #severity normal copula iid
Xfid<-rCopula(nsim, ncid) #frequency normal copula iid

Xcomons<-rCopula(nsim, normalCopula(param=1, dim =2)) #severity normal copula iid
Xcomonf<-rCopula(nsim, normalCopula(param=1, dim =2)) #severity normal copula iid
Xcounters<-rCopula(nsim, normalCopula(param=-1, dim =2)) #severity normal copula iid
Xcounterf<-rCopula(nsim, normalCopula(param=-1, dim =2)) #severity normal copula iid

# For the Model 1 in Table 7 (gamma+gamma)
#par<-result1[["optim"]][["bestmem"]] if you have run result1; see line 82
par<-c(20, 20, 986.9, 1166.8, 0.95) #otherwise, here is the result1
copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
nc <- normalCopula(param=copulaParam, dim =2)

X <- rCopula(nsim, nc) #severity normal copula 
X_md[,1] <- qgamma(X[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1 from estimated copula
X_md[,2] <- qgamma(X[,2],scale = par[3]/par[1], shape=par[1]) #severity of X2 from estimated copula
freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1 from estimated copula
freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2 from estimated copula

X_mdid[,1] <- qgamma(Xid[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1 from independent copula
X_mdid[,2] <- qgamma(Xid[,2],scale = par[3]/par[1], shape=par[1]) #severity of X2 from independent copula
freq_mdid[,1] <- qpois(Xfid[,1],lambda=0.4) #frequence of X1 illness from independent copula
freq_mdid[,2] <- qpois(Xfid[,2],lambda=0.25) #frequence of X2 accident from independent copula

SumXY[,1] <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim from estimated copula
SumXYid[,1] <- X_mdid[,1]*freq_mdid[,1] + X_mdid[,2]*freq_mdid[,2] #total claim from independent copula


# For the Model 2 in Table 7 (gamma + pareto)
#par<-result2[["optim"]][["bestmem"]]
par<-c(2.01, 20, 1049.8, 1182, 0.4)
copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
nc <- normalCopula(param=copulaParam, dim =2)
X <- rCopula(nsim, nc) #severity normal copula 
X_md[,1] <- qgamma(X[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1 from estimated copula
X_md[,2] <- qparetoII(X[,2],scale = par[3]*(par[1]-1), shape=par[1])  #severity of X2 from estimated copula
freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1 from estimated copula
freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2 from estimated copula

X_mdid[,1] <- qgamma(Xid[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1 from independent copula
X_mdid[,2] <- qparetoII(Xid[,2],scale = par[3]*(par[1]-1), shape=par[1])  #severity of X2 from independent copula
freq_mdid[,1] <- qpois(Xfid[,1],lambda=0.4) #frequence of X1 from independent copula
freq_mdid[,2] <- qpois(Xfid[,2],lambda=0.25) #frequence of X2 from independent copula
SumXY[,2] <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim from estimated copula
SumXYid[,2] <- X_mdid[,1]*freq_mdid[,1] + X_mdid[,2]*freq_mdid[,2] #total claim from independent copula

# For the Model 3 in Table 7 (gamma+weibull)
#par<-result3[["optim"]][["bestmem"]]
par<-c(10, 15, 969, 1177.3, 0.95)
copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
nc <- normalCopula(param=copulaParam, dim =2)
X <- rCopula(nsim, nc) #severity normal copula 
X_md[,1] <- qgamma(X[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1  from estimated copula
X_md[,2] <- qweibull(X[,2],scale = par[3]/gamma(1/par[1]+1), shape=par[1]) #severity of X2 from estimated copula
freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1 from estimated copula
freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2 from estimated copula

X_mdid[,1] <- qgamma(Xid[,1],scale = par[4]/par[2], shape=par[2]) #severity of X1 from independent copula
X_mdid[,2] <- qweibull(Xid[,2],scale = par[3]/gamma(1/par[1]+1), shape=par[1]) #severity of X2 from independent copula
freq_mdid[,1] <- qpois(Xfid[,1],lambda=0.4) #frequence of X1 from independent copula
freq_mdid[,2] <- qpois(Xfid[,2],lambda=0.25) #frequence of X2 from independent copula

SumXY[,3] <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim from estimated copula
SumXYid[,3] <- X_mdid[,1]*freq_mdid[,1] + X_mdid[,2]*freq_mdid[,2] #total claim from independent copula

# For the Model 4 in Table 7 (weibull+weibull)
#par<-result4[["optim"]][["bestmem"]]
par<-c(15, 15, 1070.3, 1131.7, 0.95)
copulaParam <- iRho(normalCopula(param=1, dim = 2), rho=par[5])
nc <- normalCopula(param=copulaParam, dim =2)
X <- rCopula(nsim, nc) #severity normal copula 
X_md[,1] <- qweibull(X[,1],scale = par[4]/gamma(1/par[2]+1), shape=par[2]) #severity of X1 from estimated copula
X_md[,2] <- qweibull(X[,2],scale = par[3]/gamma(1/par[1]+1), shape=par[1]) #severity of X2 from estimated copula
freq_md[,1] <- qpois(Xf[,1],lambda=0.4) #frequence of X1 from estimated copula
freq_md[,2] <- qpois(Xf[,2],lambda=0.25) #frequence of X2 from estimated copula

X_mdid[,1] <- qweibull(Xid[,1],scale = par[4]/gamma(1/par[2]+1), shape=par[2]) #severity of X1 from independent copula
X_mdid[,2] <- qweibull(Xid[,2],scale = par[3]/gamma(1/par[1]+1), shape=par[1]) #severity of X2 from independent copula
freq_mdid[,1] <- qpois(Xfid[,1],lambda=0.4) #frequence of X1 from independent copula
freq_mdid[,2] <- qpois(Xfid[,2],lambda=0.25) #frequence of X2 from independent copula

SumXY[,4] <- X_md[,1]*freq_md[,1] + X_md[,2]*freq_md[,2] #total claim from estimated copula
SumXYid[,4] <- X_mdid[,1]*freq_mdid[,1] + X_mdid[,2]*freq_mdid[,2] #total claim from independent copula

dddog<- 0:100*20 # combined deductibles
dddog[1]<-0.000001
dddog[2]<-0.0001
dddog[3]<-0.01
Ldog<-matrix(0,length(dddog),4)->Ldogid
Ldogcommon<-matrix(0,length(dddog),1)->Ldogcounter

for (i in 1:4){ # for Model 1 to Model 4
  for (j in 1:length(dddog)) {
    d = dddog[j]
    Ldog[j,i]  <- mean(pmin(SumXY[,i],d)) # deductible costs copula from model
    Ldogid[j,i]  <- mean(pmin(SumXYid[,i],d)) # deductible costs from independent model
  }
}

# Figure 7 of the paper

plot(dddog,Ldog[,1]/Ldogid[,1]-1,type="l",lty=1,col="black",xlab="Combined Deductible($)",ylab="k",main="difference in cost of deductible",ylim=c(-0.3,0))
lines(dddog,Ldog[,2]/Ldogid[,2]-1,col="red")
lines(dddog,Ldog[,3]/Ldogid[,3]-1,col="blue")
lines(dddog,Ldog[,4]/Ldogid[,4]-1,col="green")
legend("bottomright", legend=c("gamma+gamma","gamma+pareto","gamma+weibull","weibull+weibull"),col=c("black","red","blue","green"),lty = c(1,1), cex=0.85)






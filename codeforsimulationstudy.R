# This file shows how to get simulated data set named "simulationstudy.RData".
# Now we show how to simulate data for the Section 4 in the paper "Section 4 SIMULATION STUDY"

library(nlme) 
library(VGAM)
library(copula)
library(statmod)
library(numDeriv)

nsim <- 1000000 # number of simulation
set.seed(2020202)
nv <- 2 # dimension of copula is two
dparm<- 0:100/10 # range of combined deductibles 0:1000
dparm[1]<- 0.00001
dparm[2]<- 0.00005
dparm[8]<- 750/1000

rhoparam<- c(-.8,-.4,0,0.4,.8) # rhos of copula
udmax <-1000 #upper bound of the uniform distribution
shape<-matrix(c(1  ,2.5,2  ,3   ,0.5  ,1.5),nrow=3,byrow=TRUE) #shape paremeters
scale<-matrix(c(500,200,400,1200,300,400/gamma(1+1/shape[3,2])),nrow=3,byrow=TRUE) #scale paremeters
X_gamma   <- matrix(0,nsim,nv) -> X_pareto -> X_weibull
ProbSum <- matrix(0,length(dparm),length(rhoparam)) -> LimitExp -> ProbSum.2 ->LimitExpApprox.1 -> LimitExpApprox.2 ->ProbSum_g-> LimitExp_g -> ProbSum.2g ->LimitExpApprox.1g -> LimitExpApprox.2g ->ProbSum_p-> LimitExp_p -> ProbSum.2p ->LimitExpApprox.1p -> LimitExpApprox.2p -> ProbSum_w -> LimitExp_w -> ProbSum.2w ->LimitExpApprox.1w -> LimitExpApprox.2w->Probtailn->ProbtailT->ProbtailC->ProbtailF->ProbtailG
SumXY<- matrix(0,nsim,length(rhoparam))->SumXY_g->SumXY_p->SumXY_w->SumDed->SumDed_g->SumDed_p->SumDed_w
EY_g <- scale[1,1]* shape[1,1] + scale[1,2]* shape[1,2]
EY_p <- scale[2,1]/(shape[2,1]-1) + scale[2,2] /(shape[2,2]-1)
EY_w <- scale[3,1]* gamma(1+1/shape[3,1]) + scale[3,2]* gamma(1+1/shape[3,2])
covXY<-matrix(0,4,length(rhoparam))
# for different rhos (columns)

for (j in 1:length(rhoparam)) { #for five values of rho: -.8,-.4,0,0.4,.8
  # Code for different copulas
  # FrankParam <- iRho(frankCopula(param=1, dim = 2), rho=rhoparam[j])
  # nc <- frankCopula(param = FrankParam, dim = 2)
  
  # GumbelParam <- iRho(gumbelCopula(dim = 2), rho=rhoparam[j])
  # nc <- gumbelCopula(param = GumbelParam, dim = 2)  
  
  # nc <- tCopula(param=rhoparam[j], dim = 2, df = 0.5, dispstr = "ex")   
  copulaParam <- iRho(normalCopula(param=1, dim = nv), rho=rhoparam[j])
  nc <- normalCopula(param=copulaParam, dim =nv)
  X <- rCopula(nsim, nc)
  
  X_gamma[,1] <- qgamma(X[,1],scale = scale[1,1], shape=shape[1,1])
  X_gamma[,2] <- qgamma(X[,2],scale = scale[1,2], shape=shape[1,2])
  X_pareto[,1] <-  qparetoII(X[,1],scale = scale[2,1], shape=shape[2,1])
  X_pareto[,2] <-  qparetoII(X[,2],scale = scale[2,2], shape=shape[2,2])
  X_weibull[,1] <- qweibull(X[,1],scale = scale[3,1], shape=shape[3,1])
  X_weibull[,2] <- qweibull(X[,2],scale = scale[3,2], shape=shape[3,2])
  covXY[1,j]<-cov(udmax*X[,1],udmax*X[,2])
  covXY[2,j]<-cov(X_gamma[,1],X_gamma[,2])
  covXY[3,j]<-cov(X_pareto[,1],X_pareto[,2])
  covXY[4,j]<-cov(X_weibull[,1],X_weibull[,2])
  #for different combined deductible level (rows)
  for (i in 1:length(dparm)) {
    d = dparm[i]
    # Uniform(0,1000), E(X)=500
    SumXY[,j] <- udmax*X[,1] + udmax*X[,2]
    SumDed[,j] <- pmin(SumXY[,j],dparm[i]*udmax) #price of deductible
    Probtailn[i,j] <- mean(1*(X[,1] <= dparm[i] & X[,2] <= dparm[i])) # C(u,u)
    ProbSum[i,j]   <- mean(1*(SumXY[,j] <= dparm[i]*udmax)) # F(d)
    
    ProbSum.2[i,j] <- mean(1*(SumXY[,j] <= dparm[i]*udmax/2)) # F(d/2)
    LimitExp[i,j]  <- mean(pmin(SumXY[,j],dparm[i]*udmax)) # deductible costs
    LimitExpApprox.1[i,j]  <- d*udmax*(1-ProbSum[i,j])+ (d*udmax/2)*(ProbSum[i,j]) #Approximation 1
    LimitExpApprox.2[i,j]  <- d*udmax*(1-ProbSum[i,j])+ (3*d*udmax/4)*(ProbSum[i,j]-ProbSum.2[i,j])+ (d*udmax/4)*(ProbSum.2[i,j])#Approximation 2
    
    # Gamma(shape = a , scale = s); E(X) = a*s
    SumXY_g[,j] <- X_gamma[,1] + X_gamma[,2] 
    SumDed_g[,j]<- pmin(SumXY_g[,j],dparm[i]*EY_g) #price of deductible
    ProbSum_g[i,j]   <- mean(1*(SumXY_g[,j] <= dparm[i]*EY_g)) # F(d)
    ProbSum.2g[i,j] <- mean(1*(SumXY_g[,j]<= dparm[i]/2*EY_g)) # F(d/2)
    LimitExp_g[i,j]  <- mean(pmin(SumXY_g[,j],dparm[i]*EY_g)) # deductible costs
    LimitExpApprox.1g[i,j]  <- d*EY_g*(1-ProbSum_g[i,j])+ (d*EY_g/2)*(ProbSum_g[i,j]) #Approximation 1
    LimitExpApprox.2g[i,j] <- d*EY_g*(1-ProbSum_g[i,j])+ (3*d*EY_g/4)*(ProbSum_g[i,j]-ProbSum.2g[i,j])+ (d*EY_g/4)*(ProbSum.2g[i,j])#Approximation 2
    
    # Pareto(shape = a , scale = s); E(X) = s/(a-1)
    SumXY_p[,j] <- X_pareto[,1] + X_pareto[,2]
    SumDed_p[,j] <- pmin(SumXY_p[,j],dparm[i]*EY_p) #price of deductible
    ProbSum_p[i,j]   <- mean(1*(SumXY_p[,j] <= dparm[i]*EY_p)) # F(d)
    ProbSum.2p[i,j] <- mean(1*(SumXY_p[,j] <= dparm[i]/2*EY_p)) # F(d/2)
    LimitExp_p[i,j]  <- mean(pmin(SumXY_p[,j],dparm[i]*EY_p)) # deductible costs
    LimitExpApprox.1p[i,j]<- d*EY_p*(1-ProbSum_p[i,j])+ (d*EY_p/2)*(ProbSum_p[i,j]) #Approximation 1
    LimitExpApprox.2p[i,j] <- d*EY_p*(1-ProbSum_p[i,j])+ (3*d*EY_p/4)*(ProbSum_p[i,j]-ProbSum.2p[i,j])+ (d*EY_p/4)*(ProbSum.2p[i,j])#Approximation 2
    
    # Weibull(shape = a , scale = s); E(X = s*(1/a)!
    SumXY_w[,j] <- X_weibull[,1] + X_weibull[,2]
    SumDed_w[,j] <- pmin(SumXY_w[,j],dparm[i]*EY_w) #price of deductible
    ProbSum_w[i,j]   <- mean(1*(SumXY_w[,j] <= dparm[i]*EY_w)) # F(d)
    ProbSum.2w[i,j] <- mean(1*(SumXY_w[,j] <= dparm[i]/2*EY_w)) # F(d/2)
    LimitExp_w[i,j]  <- mean(pmin(SumXY_w[,j],dparm[i]*EY_w)) # deductible costs
    LimitExpApprox.1w[i,j] <- d*EY_w*(1-ProbSum_w[i,j])+ (d*EY_w/2)*(ProbSum_w[i,j]) #Approximation 1
    LimitExpApprox.2w[i,j] <- d*EY_w*(1-ProbSum_w[i,j])+ (3*d/4*EY_w)*(ProbSum_w[i,j]-ProbSum.2w[i,j])+ (d*EY_w/4)*(ProbSum.2w[i,j]) #Approximation 2
  }
}

# download P(S<d) for normal
ProbSum_n <- ProbSum
ProbSum_gn <- ProbSum_g
ProbSum_pn <- ProbSum_p
ProbSum_wn <- ProbSum_w

# Define independent case first
LimitExpiid_u <- LimitExp[,3]
LimitExpiid_g <- LimitExp_g[,3]
LimitExpiid_p <- LimitExp_p[,3]
LimitExpiid_w <- LimitExp_w[,3]

# Check the number of simulation, see Section 8 of supplementary material.
multip<-qnorm(0.975)^2*250^2
NoS=multip*t(as.matrix(c(sd(SumXY[,3])^2/mean(SumXY[,3])^2,sd(SumXY_g[,3])^2/mean(SumXY_g[,3])^2,sd(SumXY_p[,3])^2/mean(SumXY_p[,3])^2,sd(SumXY_w[,3])^2/mean(SumXY_w[,3])^2)))
(nsim>NoS)

# Then you are good to reproduce plots in the paper. For example, Figure 8.
upper3<-11
par(fig=c(0,0.52,0.65,1), new=TRUE)
plot(dparm[1:upper3]*udmax,LimitExp[1:upper3,3], ylab="Deductible Cost",xlab="",type = "l",main=" Deductible Costs",ylim=c(0,1.01*udmax))
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,4],  col="red", lty=1)
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,5],  col="red", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,1],  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,2],  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
legend("topleft", legend=c("rho=-0.8","rho=-0.4","rho=0","rho=0.4","rho=0.8"),col=c("blue","blue","black","red","red"),lty = c(2,1,1,1,2), cex=0.65,bty="n")

par(fig=c(0.48,1,0.65,1), new=TRUE)
plot(dparm[1:upper3]*udmax,LimitExp[1:upper3,5]/LimitExp[1:upper3,3]-1, lty=2, type="l", col="red", main="Cost Relativities", ylab="CR",xlab="",ylim = c(-0.15, 0.15))
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,4]/LimitExp[1:upper3,3]-1,  col="red", lty=1)
lines(dparm[1:upper3]*udmax,rep(0,length(dparm[1:upper3]*udmax)),  col="black")
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,1]/LimitExp[1:upper3,3]-1,  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp[1:upper3,2]/LimitExp[1:upper3,3]-1,  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

mtext("Marginals 1", side = 3,line =-3.5,outer = TRUE, cex = 0.75)
# gamma
par(fig=c(0,0.52,0.435,0.79), new=TRUE)
plot(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,3], type = "l",ylab="Deductible Cost",xlab="",ylim=c(0,1.01)*udmax)
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,4],  col="red", lty=1)
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,5],  col="red", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,1],  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,2],  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

par(fig=c(0.48,1,0.435,0.79), new=TRUE) 
plot(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,5]/LimitExp_g[1:upper3,3]-1, lty=2,ylab="CR",xlab="", type="l", col="red", ylim = c(-0.15, 0.15))
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,4]/LimitExp_g[1:upper3,3]-1,  col="red", lty=1)
lines(dparm[1:upper3]*udmax,rep(0,length(dparm[1:upper3]*udmax)),  col="black")
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,1]/LimitExp_g[1:upper3,3]-1,  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_g[1:upper3,2]/LimitExp_g[1:upper3,3]-1,  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
mtext("Marginals 2", side = 3,line =-13.4,outer = TRUE, cex = 0.75)
# Pareto
par(fig=c(0,0.52,0.22,0.575), new=TRUE)
plot(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,3], type = "l", ylab="Deductible Cost", xlab = "",ylim=c(0,1.01)*udmax)
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,4],  col="red", lty=1)
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,5],  col="red", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,1],  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,2],  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

par(fig=c(0.48,1,0.22,0.575), new=TRUE) 
plot(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,5]/LimitExp_p[1:upper3,3]-1, lty=2, type="l",ylab="CR", col="red", xlab = "", ylim = c(-0.15, 0.15))
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,4]/LimitExp_p[1:upper3,3]-1,  col="red", lty=1)
lines(dparm[1:upper3]*udmax,rep(0,length(dparm[1:upper3]*udmax)),  col="black")
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,1]/LimitExp_p[1:upper3,3]-1,  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_p[1:upper3,2]/LimitExp_p[1:upper3,3]-1,  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
mtext("Marginals 3", side = 3,line =-22.8,outer = TRUE, cex = 0.75)
#Weibull
par(fig=c(0,0.52,0,0.36), new=TRUE) 
plot(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,3], type = "l",ylab="Deductible Cost", xlab = "Combined Deductible: d ($)",ylim=c(0,1.01)*udmax)
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,4],  col="red", lty=1)
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,5],  col="red", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,1],  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,2],  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
par(fig=c(0.48,1,0,0.36), new=TRUE) 
plot(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,5]/LimitExp_w[1:upper3,3]-1, lty=2, type="l",ylab="CR", col="red", xlab = "Combined Deductible: d ($)", ylim = c(-0.15, 0.15))
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,4]/LimitExp_w[1:upper3,3]-1,  col="red", lty=1)
lines(dparm[1:upper3]*udmax,rep(0,length(dparm[1:upper3]*udmax)),  col="black")
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,1]/LimitExp_w[1:upper3,3]-1,  col="blue", lty=2)
lines(dparm[1:upper3]*udmax,LimitExp_w[1:upper3,2]/LimitExp_w[1:upper3,3]-1,  col="blue", lty=1)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
mtext("Marginals 4", side = 3,line =-32.6,outer = TRUE, cex = 0.75)


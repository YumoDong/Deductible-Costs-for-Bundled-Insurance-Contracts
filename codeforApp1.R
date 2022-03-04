#This file shows how to get simulated data set named "simulation1.RData".

load("pairmat.RData") # the dependence parameters for tweedies; see Table 2 of the paper.
# This is the estimated result from Jed's paper "Multivariate Frequency-Severity Regression Models in Insurance".

load("outsamplepolicyholders.RData") # this is the order of policyholders for each type coverage in the test set
# Check his personal website for the paper, data, and code. https://sites.google.com/a/wisc.edu/jed-frees/

# Now we show how to simulate data for the first application in the paper "Application 1: Wisconsin Property Fund (Commercial)"

library(copula)
B=500000   #number of simulation is 500000. Feel free to change it. But keep in mind that a lower number of simulation may give you inaccurate result in the end.
set.seed(12345)

uu <- rCopula(B,normalCopula(P2p(pairmat),dim=6,dispstr="un")) #500000 simulations from the six-dimension Gaussian copula model.

iidmatrix<-matrix(0,6,6)
diag(iidmatrix)<-1
uuiid<-rCopula(B, normalCopula(P2p(iidmatrix),dim=6,dispstr="un")) #500000 simulations from the six-dimension independent copula model.


npol<-24 # the number of policyholders we are going to randomly sample from all policyhoolders.
dparm<- 0:200*50 # the values of combined deductible from 0 to 10000.
dparm[1]<-0.000001 # revise the first three combined deductibles, this is for checking the results when combined deductible approachs 0.
dparm[2]<-0.0001
dparm[3]<-0.01
LimitExp<-matrix(0,length(dparm),npol)->LimitExpiid -> Epre -> Epreiid
simall<-matrix(0,nrow = B,ncol = npol)->simallid->simcar_dep->simcar_iid->simallcomo
sixcov<-intersect(intersect(intersect(BCyout,IMyout),intersect(COyout,CNyout)),intersect(POyout,PNyout)) # there are 244 policyholders who have all six coverages
set.seed(2020)
sample24<-sample(sixcov,npol,replace = FALSE) #randomly sample 24 policyholders from all 244 poliyholders, check Figure 1 of the paper.

#check 24 policyholders
starttime<-proc.time()
for (m in c(1:24)) {
  i=sample24[m]
  #simulation of 500000 claims for six lines from the Gaussian copula model
  sim1 <- qtweedie(uu[,1], xi=out1$xi.max, mu=fittedout[i,1], phi=out1$phi.max) #BC
  sim2 <- qtweedie(uu[,2], xi=out2$xi.max, mu=fittedout[i,2], phi=out2$phi.max) #IM
  sim3 <- qtweedie(uu[,3], xi=out3$xi.max, mu=fittedout[i,3], phi=out3$phi.max) #PN
  sim4 <- qtweedie(uu[,4], xi=out4$xi.max, mu=fittedout[i,4], phi=out4$phi.max) #PO
  sim5 <- qtweedie(uu[,5], xi=out5$xi.max, mu=fittedout[i,5], phi=out5$phi.max) #CN
  sim6 <- qtweedie(uu[,6], xi=out6$xi.max, mu=fittedout[i,6], phi=out6$phi.max) #CO
  simall[,m] <- sim1 + sim2 + sim3 + sim4 + sim5 + sim6 #total calim
  
  #simulation of 500000 claims for six lines from the independent copula model
  sim1id <- rtweedie(B, xi=out1$xi.max, mu=fittedout[i,1], phi=out1$phi.max) #BC
  sim2id <- rtweedie(B, xi=out2$xi.max, mu=fittedout[i,2], phi=out2$phi.max) #IM
  sim3id <- rtweedie(B, xi=out3$xi.max, mu=fittedout[i,3], phi=out3$phi.max) #PN
  sim4id <- rtweedie(B, xi=out4$xi.max, mu=fittedout[i,4], phi=out4$phi.max) #PO
  sim5id <- rtweedie(B, xi=out5$xi.max, mu=fittedout[i,5], phi=out5$phi.max) #CN
  sim6id <- rtweedie(B, xi=out6$xi.max, mu=fittedout[i,6], phi=out6$phi.max) #CO
  simallid[,m] <- sim1id + sim2id + sim3id + sim4id + sim5id + sim6id #total calim
}

endtime<-proc.time()
(timetaken<- -starttime+endtime) # it takes a long time since we set number of simulation.


starttime2<-proc.time()

#The deductible costs of 24 policyhoders from the copula model and the independent model under all combined deductible levels
for (m in c(1:24)) {
  for (j in 1:(length(dparm))) {
    d = dparm[j]
    LimitExp[j,m]  <- mean(pmin(simall[,m],d)) # the deductible costs of copula model
    LimitExpiid[j,m]  <- mean(pmin(simallid[,m],d)) #the deductible costs of independent model
    Epre[j,m]  <-  mean(simall[,m])- mean(pmin(simall[,m],d)) # the premiums of copula model
    Epreiid[j,m]  <- mean(simallid[,m])- mean(pmin(simallid[,m],d)) # the premiums of independent model
  }
}

stopCluster(cl)
endtime2<-proc.time()
(timetaken2<- -starttime2+endtime2)

# That's all for the simulation if you use dolloar value of combined deductible; see the left panel of Figure 3.
# If you want to use d/P% measure of deductible cost. The idea is the same. 
# You only need to divide combined deductible by the Premium_0 and repeat the above codes. The Premium_0 is stored in Epre and Epreiid.


# Evaluate Auditprob fitted model
load("fitAuditProb_Results.rdata")
require(gtools)
require(lattice)
require(coda)
require(boot)
require(rjags)
require(runjags)
require(parallel)
require(doParallel)
require(ggplot2)
library(gridExtra)
# Create expected logit_Pi and std dev
reps = dim(post)[1]
Nsmp = 1000
Detect = y
iii = sample(reps,Nsmp)
alpha = post[iii,startsWith(vn,"alpha")]
phi = post[iii,startsWith(vn,"phi[")]  
B = post[iii,startsWith(vn,"B")]  
DNNprob_reps = matrix(rep(DNNprob,Nsmp),byrow = TRUE,ncol = Nobs)
tmp = matrix(nrow = Nsmp,ncol = Nobs)
for (i in 1:Nsmp){
  tmp[i,] = DNNprob_reps[i,]*phi[i,siteN]
}
Lgt_Prob_reps = -8 + alpha[,siteN] + tmp + t(X%*%t(B))
rm(tmp)
Lgt_Prob_Mn = colMeans(Lgt_Prob_reps)
Lgt_Prob_sd = apply(Lgt_Prob_reps, 2, sd)
Expect_Prob = data.frame(Lgt_Prob_Mn=Lgt_Prob_Mn,Lgt_Prob_sd=Lgt_Prob_sd,
                         EstProb = inv.logit(Lgt_Prob_Mn),
                         site = siteN,DNNprob=DNNprob,
                         Detect=Detect)
par(mfrow=c(2,1))
ii = which(Detect==1); ii = ii[1]
plot(density(inv.logit(Lgt_Prob_reps[,ii])),main= "ExampleA: 2sec interval with Detection = 1, Estimated Prob",
     xlab = "Probability of Call Detection",ylab = "Density")
ii = which(Detect==0); ii = ii[1]
plot(density(inv.logit(Lgt_Prob_reps[,ii])),main= "ExampleB: 2sec interval with Detection = 0, Estimated Prob",
     xlab = "Probability of Call Detection",ylab = "Density")



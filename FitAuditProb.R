FitAuditProb <- function(nAuditsamps,df,Sitelist){
# nAuditsamps = c(1000) # *** MUST BE LESS THAN Nobs/Nsite, or else set to "ALL" 
# Fit Bayesian model for probability of a call from raw data using sample of audited data
# Load librar ies -------------------------------------------------------------
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
# load("Sample_raw_dat7s.rdata")
Nobs = dim(df)[1]
Nsite = dim(Sitelist)[1]
fitmodel = "MCMCprobdetect.jags" # JAGS file containing model to be solved
savename = paste0("./files/fitAuditProb_Results7_",nAuditsamps,"_AudPerSite.rdata")
savename2 = paste0("./files/fitAuditProb_Small_7s_",nAuditsamps,"_AudPerSite.rdata")
attach(df)
X = cbind(flux_sensitive,level_absolute,click,burst)
DNNprob = prob
y = Detect
siteN = site
detach(df)
Neff = dim(X)[2]
if (nAuditsamps!="ALL"){
  iii = numeric()
  for (i in 1:Nsite){
    iii = c(iii,sample(which(df$site==i),nAuditsamps,replace = FALSE))
  }
  X = X[iii,]
  DNNprob = DNNprob[iii]
  y = y[iii]
  siteN = siteN[iii]
  Nobs1 = length(iii)
}else{
  Nobs1 = Nobs
}
lgDNN = pmax(-15,logit(DNNprob))
hist(lgDNN)
# Priors and initial values (to speed up fitting - can be set more vague or no inits)
# LowB = c(4, .001,-1.5,-3,1)
# HighB= c(6, .005,-0.5,-1,3)
sigPri = c(1,2.5,2.5,2.5)
tauPri = 1/sigPri^2
# Note: if 5 variables in Matrix X, and 5 "B" params, then calc eqn using inner product
#  Prob = inv.logit(-8 + alpha[i] + phi[i]*DNNprob + X%*%B) 
# where eps = either actual random effect for site i, OR random effect from sigma
# Set up Jags inputs -----------------------------------------------------
#
jags.data <- list(y=y,Nobs=Nobs1,Nsite=Nsite,site=siteN,X=X,DNNprob=DNNprob,Neff=Neff,tauB=tauPri) # 
#
inits <- function() list(sigA = runif(1,1.5,2.5),sigH = runif(1,.2,.4),
                         phi0 = runif(1,4,6), alph0 = runif(1,-6,-5)) #,B = runif(Neff, LowB, HighB) B0 = runif(1,-11,-9), 
#
params <- c("sigA","sigH","B","phi0","alph0","phi","alpha") # N Dispers

nsamples <- 500
nt <- 1
nb <- 3000
cores = detectCores()
ncore = min(20,cores-1)
nc <- ncore
cl <- makeCluster(ncore)
#
# Run JAGS to fit model---------------------------------------------
#
out <- run.jags(data = jags.data, 
                inits=inits,
                monitor = params, 
                model = fitmodel, 
                n.chains = nc, 
                thin = nt, 
                sample = nsamples, 
                burnin = nb,
                method="rjparallel", cl=cl)
#
# Diagnostic plots -------------------------------------------------
#
stopCluster(cl)
post = rbind(out$mcmc[[1]], out$mcmc[[2]])
for (i in 3:nc){
  post = rbind(post, out$mcmc[[i]])
}
sumstats = summary(out)
vn = row.names(sumstats)
#
plot(out, vars = c("sigA","sigH","phi0","alph0","B","alpha","phi"), 
     plot.type = c("trace", "histogram"), layout = c(1,2))
# ----
reps = dim(post)[1]
Nsmp = 5000
attach(df)
X = cbind(flux_sensitive,level_absolute,click,burst)
DNNprob = prob
y = Detect
siteN = site
detach(df)
Detect = y
iii = sample(reps,Nsmp)
# Genertate expected Logit_P values and uncertainty (means and sd) from posterior
#  THIS CODE OR SOMETHING LIKE IT CAN BE APPLIED TO ALL 2S WINDOWS:
alpha = post[iii,startsWith(vn,"alpha")]
phi = post[iii,startsWith(vn,"phi[")]  
B = post[iii,startsWith(vn,"B")]  
DNNprob_reps = matrix(rep(DNNprob,Nsmp),byrow = TRUE,ncol = Nobs)
tmp = matrix(nrow = Nsmp,ncol = Nobs)
for (i in 1:Nsmp){
  tmp[i,] = DNNprob_reps[i,]*phi[i,siteN]
}
Lgt_Prob_reps = alpha[,siteN] + tmp + t(X%*%t(B))
rm(tmp)
Lgt_Prob_Mn = colMeans(Lgt_Prob_reps)
Lgt_Prob_sd = apply(Lgt_Prob_reps, 2, sd); 
# Convert from Logits to probabilities with variances using delta method 
EstProb_sd = ((exp(Lgt_Prob_Mn)/(1+exp(Lgt_Prob_Mn))^2)^2)*Lgt_Prob_sd^2 
EstProb = inv.logit(Lgt_Prob_Mn)-0.5*( exp(Lgt_Prob_Mn)*(exp(Lgt_Prob_Mn)-1)/(1+exp(Lgt_Prob_Mn))^3)*Lgt_Prob_sd^2
Expect_Prob = data.frame(Lgt_Prob_Mn=Lgt_Prob_Mn,Lgt_Prob_sd=Lgt_Prob_sd,
                         EstProb = EstProb, EstProb_sd = EstProb_sd,
                         site = siteN,DNNprob=DNNprob,
                         Detect=Detect)
par(mfrow=c(2,1))
ii = which(Detect==1); ii = ii[1]
plot(density(inv.logit(Lgt_Prob_reps[,ii])),main= "ExampleA: 2sec interval with Detection = 1, Estimated Prob",
     xlab = "Probability of Call Detection",ylab = "Density")
ii = which(Detect==0); ii = ii[1]
plot(density(inv.logit(Lgt_Prob_reps[,ii])),main= "ExampleB: 2sec interval with Detection = 0, Estimated Prob",
     xlab = "Probability of Call Detection",ylab = "Density")
par(mfrow=c(1,1))
#
save(df,Sitelist,Expect_Prob,out,post,sumstats,nAuditsamps,file=savename) 
save(df,Sitelist,Expect_Prob,sumstats,nAuditsamps,file=savename2) 
#
result <- list(NAudits=nAuditsamps,
               sumstats = sumstats,
               savename = savename,
               savenameSH = savename2)
return(result)
}
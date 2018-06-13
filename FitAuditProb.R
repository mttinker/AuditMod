# Fit Bayesian model for probability of a call from raw data using sample of audited data
# Load libraries -------------------------------------------------------------
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
load("Sample_raw_dat7s.rdata")
nAuditsamps = c(1000) # *** MUST BE LESS THAN Nobs/Nsite, or else set to "ALL" 
fitmodel = "MCMCprobdetect.jags" # JAGS file containing model to be solved
savename = paste0("fitAuditProb_Results_",nAuditsamps,"_AudPerSite.rdata")
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
inits <- function() list(sigS = runif(1,1.5,2.5),sigP = runif(1,.2,.4),
                         phi0 = runif(1,4,6)) #,B = runif(Neff, LowB, HighB) B0 = runif(1,-11,-9), 
#
params <- c("sigS","sigP","B","phi0","phi","alpha") # N Dispers

nsamples <- 500
nt <- 1
nb <- 2000
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
plot(out, vars = c("sigS","sigP","phi0","B","alpha","phi"), plot.type = c("trace", "histogram"), layout = c(1,2))
# ----
reps = dim(post)[1]
Nsmp = 1000
attach(df)
X = cbind(flux_sensitive,level_absolute,click,burst)
DNNprob = prob
y = Detect
siteN = site
detach(df)
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
save.image(file=savename)
par(mfrow=c(1,1))
## Estmate error for Audited detection vs estimate prob by site -------------------
df1 = df
sampTS = 5
nsamp = 100
sampH = round(.33*nsamp)
sampL = nsamp-sampH
reps = 5000
Naudit = round(Nobs1/Nsite)
site = numeric(length = Nobs)
par(mfrow=c(2,2))
for (s in 1:Nsite){
  iii = which(df1$location==Sitelist[s,1])
  site[iii] = s
  df = df1[iii,]
  df$lgtprobfit = Expect_Prob$Lgt_Prob_Mn[iii]
  df$lgtprobse = Expect_Prob$Lgt_Prob_sd[iii]
  df$probfit = Expect_Prob$EstProb[iii]
  RecrdsH = which(df$prob>=0.95); ix = which((RecrdsH+sampTS)<length(iii)); RecrdsH=RecrdsH[ix]
  RecrdsL = which(df$prob<0.95); ix = which((RecrdsL+sampTS)<length(iii)); RecrdsL=RecrdsL[ix]
  SimsampT1 = matrix(data=0,nrow = sampH, ncol = sampTS)
  SimsampE1 = matrix(data=0,nrow = sampH, ncol = sampTS) 
  SimsampT2 = matrix(data=0,nrow = sampL, ncol = sampTS)
  SimsampE2  = matrix(data=0,nrow = sampL, ncol = sampTS)
  EstBias = numeric(length = reps)
  EstErr = numeric(length = reps)
  CallsT= numeric(length = reps)
  CallsE= numeric(length = reps)
  for (i in 1:reps){
    ii = sample(RecrdsH,sampH,replace = TRUE)
    for (j in 1:sampTS){
      SimsampT1[,j] = df$Detect[ii+j-1]
      SimsampE1[,j] = inv.logit(rnorm(sampH,df$lgtprobfit[ii+j-1],df$lgtprobse[ii+j-1]))
    }
    ii = sample(RecrdsL,sampL,replace = TRUE)
    for (j in 1:sampTS){
      SimsampT2[,j] = df$Detect[ii+j-1]
      SimsampE2[,j] = inv.logit(rnorm(sampL,df$lgtprobfit[ii+j-1],df$lgtprobse[ii+j-1]))
    }  
    SimsampT = rbind(SimsampT1,SimsampT2)
    SimsampE = rbind(SimsampE1,SimsampE2)
    # tmp1 = fitdist(rowSums(SimsampT),"nbinom")
    # tmp2 = fitdist(rowSums(SimsampE),"lnorm")
    # CallsT[i] = tmp1$estimate[2] 
    # CallsE[i] = exp(tmp2$estimate[1] + tmp2$estimate[2]^2/2)
    CallsT[i] = mean(rowSums(SimsampT))/sampTS*30
    CallsE[i] = mean(rowSums(SimsampE))/sampTS*30
    EstBias[i] = CallsE[i] - CallsT[i]
    EstErr[i] = (CallsE[i] - CallsT[i])^2
  }
  MnEstBias = mean(EstBias)
  EstStdErr = sqrt(sum(EstErr)/(reps-1)) 
  tmp = density(CallsT); maxy = max(tmp$y)
  plot(density(c(CallsT,CallsE)),type="n", ylim=c(0,1.5*maxy),
       main = paste0("Site ",Sitelist[s,1], "N Audits = ",format(Naudit,digits = 0),
                     ", StdErr = ", format(EstStdErr,digits = 3),
                     ", Mean Bias = ", format(MnEstBias,digits = 3)),
       xlab="Estimated Mean Call Rate", ylab="Probability Density")
  lines(density(CallsT),col="blue")
  lines(density(CallsE),col="red")
  legend(mean(CallsT), maxy, legend=c("Estimate from Audited Detections",
                                      "Probability-based Estimates"),
         col=c("blue", "red"), lty=1, cex=0.8)
}
df = df1

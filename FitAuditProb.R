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
fitmodel = "MCMCprobdetect.jags" # JAGS file containing model to be solved
load("Sample_raw_dat.rdata")
attach(df)
X = cbind(flux_sensitive,level_absolute,click,burst)
DNNprob = prob
y = Detect
siteN = site
detach(df)
Neff = dim(X)[2]
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
jags.data <- list(y=y,Nobs=Nobs,Nsite=Nsite,site=siteN,X=X,DNNprob=DNNprob,Neff=Neff,tauB=tauPri) # 
#
inits <- function() list(sigS = runif(1,1.5,2.5),sigP = runif(1,.2,.4),
                         phi0 = runif(1,4,6)) #,B = runif(Neff, LowB, HighB) B0 = runif(1,-11,-9), 
#
params <- c("sigS","sigP","B","phi0","phi","alpha") # N Dispers

nsamples <- 500
nt <- 1
nb <- 2500
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
save.image(file="fitAuditProb_Results.rdata")
# ----


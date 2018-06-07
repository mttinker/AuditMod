# Analyze effect of number of audited records on estimation precision
# 
# User inputs -----------------------------------------------------------
df = read.csv("OurData_Subset_Audited.csv")
TargetCall = c("HAPE_h")
RcrdsperTS = 30
nsamp = 5000 # Number of 15min time step records per site, hypothetical
reps = 1000  # Number of replications for simulations to estimate std errors
Nobs = length(df$rating)
AuditNs = c(250,seq(500,5000,by=500),Nobs) #Vector of possible Audit sample sizes
# 
# Load Libraries --------------------------------------------------------
#
library(gtools)
require(parallel)
require(doParallel)
require(foreach)
library(ggplot2)
library(gridExtra)
# library(fitdistrplus)
#
# Process data ----------------------------------------------------------
ratings = levels(df$rating)
ratetrue = grep(TargetCall,ratings)
ii = numeric()
for (r in 1:length(ratetrue)){
  iii = which(df$rating==ratings[ratetrue[r]]) 
  ii = c(ii,iii)
}
Detect = numeric(length = Nobs)
Detect[ii] = 1
df$Detect =Detect 
df$flux_sens2 = df$flux_sensitive^2
df$flux_lev_intxn = df$flux_sensitive*df$level_absolute
RecrdsH = which(df$prob>=0.95)
RecrdsL = which(df$prob<0.95)
df1 = df
rm(df)
Nss = length(AuditNs)
sampH = round(RcrdsperTS*(1/3)) # Proportion of high-probability samples
sampL = RcrdsperTS-sampH        # Proportion of low-probability samples
MnEstBias = numeric(length = Nss)
EstStdErr = numeric(length = Nss)
SimsampT1 = matrix(data=0,nrow = nsamp, ncol = sampH)
SimsampE1 = matrix(data=0,nrow = nsamp, ncol = sampH) 
SimsampT2 = matrix(data=0,nrow = nsamp, ncol = sampL)
SimsampE2 = matrix(data=0,nrow = nsamp, ncol = sampL)
EstBias = numeric(length = reps)
EstErr = numeric(length = reps)
CallsT = numeric(length = reps)
CallsE = numeric(length = reps)
## Estmate error for Audited detection vs estimate prob -------------------
#
Availcores = detectCores()
ncores = min(Nss,Availcores)
cl <- makeCluster(ncores)
registerDoParallel(cl)
nodeNames <- foreach(i = 1:length(cl), .combine=c) %dopar% {
  Sys.getpid()
}
# for (s in 1:Nss){
results <- foreach(q = 1:ncores) %dopar% {
  s = which(Sys.getpid()==nodeNames)
  df = df1
  Naudit = AuditNs[s]
  iii = sample(Nobs,Naudit,replace = FALSE)
  # Fit logistic model of realized probability vs predicted prob given machine learning prediction
  # (with signal quality metrics as co-variates)
  mod1 <- glm(Detect ~ prob + flux + flux_sensitive + flux_sens2 + level + level_absolute + click + burst, 
              data=df[iii,], family=binomial(link="logit"))
  # summary(mod1)
  preddat <- predict(mod1, newdata=df, se.fit=TRUE)
  df$lgtprobfit = preddat$fit
  df$lgtprobse = preddat$se.fit
  # df$probfit = inv.logit(preddat$fit)
  df$probfit = exp(preddat$fit)/(1+exp(preddat$fit))
  for (i in 1:reps){
    for (j in 1:sampH){
      ii = sample(RecrdsH,nsamp,replace = TRUE)
      SimsampT1[,j] = df$Detect[ii]
      SimsampE1[,j] = df$probfit[ii]
    }
    for (j in 1:sampL){
      ii = sample(RecrdsL,nsamp,replace = TRUE)
      SimsampT2[,j] = df$Detect[ii]
      SimsampE2[,j] = df$probfit[ii]
    }  
    SimsampT = cbind(SimsampT1,SimsampT2)
    SimsampE = cbind(SimsampE1,SimsampE2)
    CallsT[i] = mean(rowSums(SimsampT))
    CallsE[i] = mean(rowSums(SimsampE))
    EstBias[i] = CallsE[i] - CallsT[i]
    EstErr[i] = (CallsE[i] - CallsT[i])^2
  }
  MnEstBias = mean(EstBias)
  EstStdErr = sqrt(sum(EstErr)/(reps-1)) 
  LgtProb_se = mean(df$lgtprobse)
  Simresults = list(MnEstBias=MnEstBias,EstStdErr=EstStdErr,Naudit=Naudit,
                        LgtProb_se=LgtProb_se,CallsT=CallsT,CallsE=CallsE)
  return(Simresults)  
}
stopCluster(cl)
# Extract and review results --------------------------------------------
MnEstBias = numeric(length = Nss)
EstStdErr = numeric(length = Nss)
LgtProb_se = numeric(length = Nss)
CallsT = matrix(nrow = reps,ncol = Nss)
CallsE = matrix(nrow = reps,ncol = Nss)
for (s in 1:Nss){
  MnEstBias[s] = results[[s]]$MnEstBias
  EstStdErr[s] = results[[s]]$EstStdErr
  LgtProb_se[s] = results[[s]]$LgtProb_se
  CallsT[,s] = results[[s]]$CallsT
  CallsE[,s] = results[[s]]$CallsE
}
#
Auditdf = data.frame(N_Audit_Records=AuditNs,MnEstBias=MnEstBias,
                     EstPrbSE=EstStdErr,LgtProbSE=LgtProb_se)
mod1 = lm(log(LgtProbSE)~log(N_Audit_Records),data=Auditdf)
prd1 = predict(mod1,Auditdf,se.fit = TRUE,
              interval = c("confidence"),level = 0.95)
Auditdf$LgtProbSE_sm = exp(prd1$fit[,1])
Auditdf$LgtProbSE_smLO = exp(prd1$fit[,2])
Auditdf$LgtProbSE_smHI = exp(prd1$fit[,3])
mod2 = lm(log(EstPrbSE)~log(N_Audit_Records),data=Auditdf)
prd2 = predict(mod2,Auditdf,se.fit = TRUE,
               interval = c("confidence"),level = 0.95)
Auditdf$EstPrbSE_sm = exp(prd2$fit[,1])
Auditdf$EstPrbSE_smLO = exp(prd2$fit[,2])
Auditdf$EstPrbSE_smHI = exp(prd2$fit[,3])
#
plt1 <- ggplot(Auditdf, aes(x = N_Audit_Records, y = LgtProbSE)) + geom_point() +
  geom_ribbon(aes(ymin=LgtProbSE_smLO,ymax=LgtProbSE_smHI),alpha=0.3)+
  geom_line(aes(y = LgtProbSE_sm),size=1)+
  xlab("Number of Audited Records (used to build probability model)") +
  ylab("Mean SE of Logit Probabilities") +
  ggtitle("Logit Probability Estimate SE vs Audit Sample Size")
# print(plt1)
plt2 <- ggplot(Auditdf, aes(x = N_Audit_Records, y = EstPrbSE)) + geom_point() +
  geom_ribbon(aes(ymin=EstPrbSE_smLO,ymax=EstPrbSE_smHI),alpha=0.3)+
  geom_line(aes(y = EstPrbSE_sm),size=1)+
  xlab("Number of Audited Records (used to build probability model)") +
  ylab("SE of Model-Estimated Call Rate") +
  ggtitle("Call Rate Estimation SE vs Audit Sample Size")
# print(plt2)
grid.arrange(plt1, plt2, nrow=2)
#
s = 3
p1t3 = ggplot(data.frame(CallsT=CallsT[,1],CallsE=CallsE[,s])) + 
  geom_density(aes(CallsT, colour = "CallsT"),size=1) + geom_density(aes(CallsE, colour = "CallsE"),size=1) + 
  xlab("Estimated Calls per Minute") + ylab("Distribution Density") +
  ggtitle(paste0("Call Rate, Audited Detections vs Probability Estimates, N_Audit = ", AuditNs[s]))+
  scale_colour_discrete(name="Method of Estimation",
                        breaks=c("CallsT", "CallsE"),
                        labels=c("Audited Detections Only", "Model-Estimated Probabilities"))
s = 12
p1t4 = ggplot(data.frame(CallsT=CallsT[,1],CallsE=CallsE[,s])) + 
  geom_density(aes(CallsT, colour = "CallsT"),size=1) + geom_density(aes(CallsE, colour = "CallsE"),size=1) + 
  xlab("Estimated Calls per Minute") + ylab("Distribution Density") +
  ggtitle(paste0("Call Rate, Audited Detections vs Probability Estimates, N_Audit = ", AuditNs[s]))+
  scale_colour_discrete(name="Method of Estimation",
                        breaks=c("CallsT", "CallsE"),
                        labels=c("Audited Detections Only", "Model-Estimated Probabilities"))
grid.arrange(p1t3, p1t4, nrow=2)

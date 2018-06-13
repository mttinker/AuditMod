# Data Review
library(gtools)
library(fitdistrplus)
df = read.csv("OurData_Subset_Audited.csv")
TargetCall = c("HAPE_h")
# 
Nobs = length(df$rating)
ratings = levels(df$rating)
ratetrue = grep(TargetCall,ratings)
ii = numeric()
for (r in 1:length(ratetrue)){
  iii = which(df$rating==ratings[ratetrue[r]]) 
  ii = c(ii,iii)
}
Detect = numeric(length = Nobs)
Detect[ii] = 1
Prob1 = df$prob
plot(density(Prob1[Detect==1],adjust = .2),col="blue",
     main = "Estimated Probs for Positive vs Negative Detections",
     xlab = "Modeled Probability (Machine Learning)")
lines(density(Prob1[Detect==0],adjust = .2),col="red")
tmp = density(Prob1[Detect==1],adjust = .2); maxy = max(tmp$y)
legend(.1, .8*maxy, legend=c("Postitive Detection", "Negative Detection"),
      col=c("blue", "red"), lty=1, cex=0.8)
Crit = 0.95
df$Outcm = numeric(length=Nobs)
df$Outcm[which(Detect==1 & Prob1>=Crit)] = 1
df$Outcm[which(Detect==1 & Prob1<Crit)] = 2
df$Outcm[which(Detect==0 & Prob1<Crit)] = 3
df$Outcm[which(!is.na(df$rating) & Detect==0 & Prob1>=Crit)] = 4
Outcomes = c("True pos","False neg","True neg","False pos")
df$Detect =Detect 
# df$Outcm = factor(df$Outcm,labels = Outcomes)
df$flux_sens2 = df$flux_sensitive^2
df$flux_lev_intxn = df$flux_sensitive*df$level_absolute
df1 = df
Locs = levels(df$location)
NLocs = length(Locs); Nsite = NLocs
par(mfrow=c(2,2))
site = numeric(length = Nobs)
for (i in 1:NLocs){
  site[which(df1$location==Locs[i])] = i
  df = df1[which(df1$location==Locs[i]),]
  tmp = hist(df$Outcm,seq(.5,4.5,by=1),plot = FALSE)
  barplot(tmp$counts, xaxt = "n", xlab='Outcomes',
        main = paste0("Location = ",Locs[i]))
  axis(1, at=seq(.75,4.25,length.out = 4), labels=Outcomes)
}
par(mfrow=c(1,1))
df = df1
df$site = site
Sitelist = data.frame(Site = Locs, SiteN = seq(1,Nsite))
# Fit model of actual probability vs predicted prob given machine learning prediction
rm(site,Detect)
save(df,Sitelist,Outcomes,ratings,Nsite,Nobs,file="Sample_raw_dat.rdata")

# Model for ALL sites --------------------------------------------------------------
mod1 <- glm(Detect ~ prob + flux_sensitive + level_absolute + click + burst , 
            data=df, family=binomial(link="logit"))
summary(mod1)
preddat <- predict(mod1, newdata=df, se.fit=TRUE)
df$lgtprobfit = preddat$fit
df$lgtprobse = preddat$se.fit
df$probfit = inv.logit(preddat$fit)

# Plot function with mean qulity metrics
newdat = data.frame(prob = seq(0,1,by=.01),flux = rep(mean(df$flux),101),
                    flux_sensitive=rep(mean(df$flux_sensitive),101),
                    flux_sens2 = rep(mean(df$flux_sensitive)^2,101),
                    level = rep(mean(df$level),101), 
                    level_absolute = rep(mean(df$level_absolute),101),
                    click = rep(mean(df$click),101),
                    burst = rep(mean(df$burst),101))
plotdat <- predict(mod1, newdata=newdat, se.fit=TRUE)
with(df, plot(prob, Detect, type="n", 
                 ylim=c(0, 1), xlab="Predicted Prob", ylab="Prob Actual Detection"))
with(plotdat, lines(seq(0,1,by=.01), inv.logit(fit), col="blue"))
with(plotdat, lines(seq(0,1,by=.01), inv.logit(fit+1.96*se.fit), lty=2))
with(plotdat, lines(seq(0,1,by=.01), inv.logit(fit-1.96*se.fit), lty=2))
#
## Estmate error for Audited detection vs estimate prob by site -------------------
sampH = 10
sampL = 30-sampH
nsamp = 5000
reps = 500
Naudit = 1500
site = numeric(length = Nobs)
par(mfrow=c(2,2))
for (s in 1:Nsite){
  site[which(df1$location==Locs[s])] = s
  df = df1[which(df1$location==Locs[s]),]
  mod <- glm(Detect ~ prob + flux_sensitive + level_absolute + click + burst , 
             data=df, family=binomial(link="logit"))
  txt = summary(mod)
  print(txt)
  preddat <- predict(mod, newdata=df, se.fit=TRUE)
  df$lgtprobfit = preddat$fit
  df$lgtprobse = preddat$se.fit
  df$probfit = inv.logit(preddat$fit)
  RecrdsH = which(df$prob>=0.95)
  RecrdsL = which(df$prob<0.95)
  SimsampT1 = matrix(data=0,nrow = nsamp, ncol = sampH)
  SimsampE1 = matrix(data=0,nrow = nsamp, ncol = sampH) 
  SimsampT2 = matrix(data=0,nrow = nsamp, ncol = sampL)
  SimsampE2  = matrix(data=0,nrow = nsamp, ncol = sampL)
  EstBias = numeric(length = reps)
  EstErr = numeric(length = reps)
  CallsT= numeric(length = reps)
  CallsE= numeric(length = reps)
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
    # tmp1 = fitdist(rowSums(SimsampT),"nbinom")
    # tmp2 = fitdist(rowSums(SimsampE),"lnorm")
    # CallsT[i] = tmp1$estimate[2] 
    # CallsE[i] = exp(tmp2$estimate[1] + tmp2$estimate[2]^2/2)
    CallsT[i] = mean(rowSums(SimsampT))
    CallsE[i] = mean(rowSums(SimsampE))
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



# Evaluate results of different sample sizes of audited samples (power analysis)
library(ggplot2)
library(reshape2)
Text1 = c("./files/fitAuditProb_Small_7s_")
Text2 = c("_AudPerSite.rdata")
AuditSzs = c(100,300,500,750,1000)
NSS = length(AuditSzs)
FileList = character(length = NSS)
for (f in 1:NSS){
  FileList[f] = paste0(Text1,AuditSzs[f],Text2)
}
load(FileList[1])
Nsite = dim(Sitelist)[1]
MnEstDev = matrix(nrow = NSS,ncol = Nsite)
EstStdErr = matrix(nrow = NSS,ncol = Nsite)
for (f in 1:NSS){
  load(FileList[f])
  df1 = df
  sampTS = 5
  nsamp = 100
  sampH = round(.25*nsamp)
  sampL = nsamp-sampH
  reps = 5000
  Naudit = AuditSzs[f]
  site = numeric(length = dim(df)[1])
  par(mfrow=c(2,2))
  for (s in 1:Nsite){
    iii = which(df1$location==Sitelist[s,1])
    site[iii] = s
    df = df1[iii,]
    df$lgtprobfit = Expect_Prob$Lgt_Prob_Mn[iii]
    df$lgtprobse = Expect_Prob$Lgt_Prob_sd[iii]
    df$probfit = Expect_Prob$EstProb[iii]
    RecrdsH = which(df$prob>=0.5); ix = which((RecrdsH+sampTS)<length(iii)); RecrdsH=RecrdsH[ix]
    RecrdsL = which(df$prob<0.5); ix = which((RecrdsL+sampTS)<length(iii)); RecrdsL=RecrdsL[ix]
    SimsampT1 = matrix(data=0,nrow = sampH, ncol = sampTS)
    SimsampE1 = matrix(data=0,nrow = sampH, ncol = sampTS) 
    SimsampT2 = matrix(data=0,nrow = sampL, ncol = sampTS)
    SimsampE2  = matrix(data=0,nrow = sampL, ncol = sampTS)
    EstDev = numeric(length = reps)
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
      EstDev[i] = CallsE[i] - CallsT[i]
      EstErr[i] = (CallsE[i] - CallsT[i])^2
    }
    MnEstDev[f,s] = abs(mean(EstDev))
    EstStdErr[f,s] = sqrt(sum(EstErr)/(reps-1)) 
    tmp = density(CallsT); maxy = max(tmp$y)
    plot(density(c(CallsT,CallsE)),type="n", ylim=c(0,1.5*maxy),
         main = paste0("Site ",Sitelist[s,1], "N Audits = ",format(Naudit,digits = 0),
                       ", StdErr = ", format(EstStdErr[f,s],digits = 3),
                       ", Mean Dev = ", format(MnEstDev[f,s],digits = 3)),
         xlab="Estimated Mean Call Rate", ylab="Probability Density")
    lines(density(CallsT),col="blue")
    lines(density(CallsE),col="red")
    legend(mean(CallsT), maxy, legend=c("Estimate from Audited Detections",
                                        "Probability-based Estimates"),
           col=c("blue", "red"), lty=1, cex=0.8)
  }
}
par(mfrow=c(1,1))

SiteDev = as.data.frame(cbind(AuditSzs,MnEstDev)); colnames(SiteDev) = c("N_Audits",as.character(Sitelist$Site))
SiteErr = as.data.frame(cbind(AuditSzs,EstStdErr)); colnames(SiteErr) = c("N_Audits",as.character(Sitelist$Site))

MeanDev = rowMeans(MnEstDev)
MeanErr = rowMeans(EstStdErr)
dfDevMn = data.frame(x = AuditSzs, y = MeanDev)
dfErrMm = data.frame(x = AuditSzs, y = MeanErr)
lgY = log(MeanDev); lgX = log(AuditSzs); 
fit1 = lm(lgY ~ lgX )
# summary(fitExt)
newdat = data.frame(lgX = log(seq(100,1000)))
Prd = predict(fit1,newdat,se.fit = TRUE,
              interval = c("confidence"),level = 0.95)
dfplt1 = data.frame(SampSize = seq(100,1000),
                    MeanDev = exp(Prd$fit[,1]),
                    CILO = exp(Prd$fit[,2]),CIHI = exp(Prd$fit[,3]))
plt1 = (ggplot(dfplt1, aes(x= SampSize))+
             geom_ribbon(aes(ymin=CILO,ymax=CIHI),alpha=0.3)+
             geom_line(aes(y = MeanDev),size=1)+
             geom_point(data = dfDevMn, aes(x=x,y=y),size=2) +
             xlab("Number of Audited Samples") +
             ylab("Mean Deviation, abs(Estimate_Audit - Estimate_Prob)") +
             ggtitle("Deviation in Call Rate Estimates, Audited Detections vs Predicted Probs"))
print(plt1)

dfErr = data.frame(x = AuditSzs, y = MeanErr)
lgY = log(MeanErr); lgX = log(AuditSzs); 
fit2 = lm(lgY ~ lgX )
# summary(fitExt)
newdat = data.frame(lgX = log(seq(100,1000)))
Prd = predict(fit2,newdat,se.fit = TRUE,
              interval = c("confidence"),level = 0.95)
dfplt2 = data.frame(SampSize = seq(100,1000),
                    MeanError = exp(Prd$fit[,1]),
                    CILO = exp(Prd$fit[,2]),CIHI = exp(Prd$fit[,3]))
plt2 = (ggplot(dfplt2, aes(x= SampSize))+
          geom_ribbon(aes(ymin=CILO,ymax=CIHI),alpha=0.3)+
          geom_line(aes(y = MeanError),size=1)+
          geom_point(data = dfErr, aes(x=x,y=y)) +
          xlab("Number of Audited Samples") +
          ylab("Mean Standard Error") +
          ggtitle("Std Error in Call Rate Estimates, Audited Detections vs Predicted Probs"))
print(plt2)

dfplt3 <- melt(SiteDev,  id.vars = 'N_Audits', variable.name = 'Site')
plt3 = ggplot(dfplt3, aes(N_Audits,value)) + geom_line(aes(colour = Site)) +
  xlab("Number of Audited Samples") +
  ylab("Mean Deviation, abs(Estimate_Audit - Estimate_Prob)") +
  ggtitle("Deviation in Call Rate Estimates by Site, Audited Detections vs Predicted Probs")
print(plt3)

dfplt4 <- melt(SiteErr,  id.vars = 'N_Audits', variable.name = 'Site')
plt4 = ggplot(dfplt4, aes(N_Audits,value)) + geom_line(aes(colour = Site)) +
  xlab("Number of Audited Samples") +
  ylab("Mean Standard Error") +
  ggtitle("Std Error in Call Rate Estimates by Site, Audited Detections vs Predicted Probs")
print(plt4)


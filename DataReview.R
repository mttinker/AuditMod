# Data Review
library(gtools)
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
df2 = df[df$Outcm>0,]
Prob1 = Prob1[df$Outcm>0]
Detect = Detect[df$Outcm>0]
df2$Detect =Detect 
tmp = hist(df$Outcm,seq(.5,4.5,by=1),plot = FALSE)
df2$Outcm = factor(df2$Outcm,labels = Outcomes)
df1 = df
df=df2
barplot(tmp$counts, xaxt = "n", xlab='Outcomes',
        main = paste0("Critical value = ",Crit ))
axis(1, at=seq(.75,4.25,length.out = 4), labels=Outcomes)
# Fit model of actual probability vs predicted prob given machine learning prediction
mod1 <- glm(Detect ~ prob + flux + flux_sensitive + level + level_absolute + click + burst, 
            data=df, family=binomial(link="logit"))
summary(mod1)
preddat <- predict(mod1, newdata=df, se.fit=TRUE)
df$lgtprobfit = preddat$fit
df$lgtprobse = preddat$se.fit
df$probfit = inv.logit(preddat$fit)

# Plot function with mean qulity metrics
newdat = data.frame(prob = seq(0,1,by=.01),flux = rep(mean(df$flux),101),
                    flux_sensitive=rep(mean(df$flux_sensitive),101),
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

# Estimate error based on 500 random samples of 30 iterated 100 times
# (with estimated prob vs audited detections)
samp = 30
reps1 = 500
reps2 = 100
trueM = numeric(length = reps1*reps2)
estM = numeric(length = reps1*reps2)
dev2 = numeric(length = reps1*reps2)
r = 0
for (r1 in 1:reps1){
  ii = sample(Nobs,samp,replace = TRUE)
  for (r2 in 1:reps2){
    r = r+1
    trueM[r] = sum(df$Detect[ii])
    estM[r] = sum(rbinom(samp,1,inv.logit(rnorm(samp,df$lgtprobfit[ii],df$lgtprobse[ii]))))
    dev2[r] = (estM[r]-trueM[r])^2
  }
}
dev = sqrt(mean(dev2))
plot(density(trueM),main = paste0("StdErr = ",dev),xlab="Mean Detection")
lines(density(estM))

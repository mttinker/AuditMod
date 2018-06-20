library(readxl)
library(lubridate)
library(plyr)
library(stringr)
library(gtools)
library(coda)
library(lattice)
library(rjags)
library(coda)
library(mcmcplots)
library(jagsUI)
library(parallel)
library(doParallel)
library(fitdistrplus)
# LOAD RAW DATA --------------------------------------------------------------------
loadfile1 = c("./files/fitAuditProb_Small_7s_1000_AudPerSite.rdata")
load(loadfile1)
df1 = df
# Process Data files ---------------------------------------------------------------
iii = which(df$TimeStepN<=50)
df = df1[iii,]
Night1 = min(df$NightOfYear)
WeekN =ceiling((df$NightOfYear-Night1+1)/7)
DayN = (df$NightOfYear-Night1+1)
TS = df$TimeStepN
DNN = df$prob
Flux = df$flux_sensitive
LvAb = df$level_absolute
Clck = df$click
Brst = df$burst
NWeeks = max(WeekN)
NSite = dim(Sitelist)[1]
Illu = df$Moonup*df$Illu
Moon = 2*Illu-1

Prob = Expect_Prob$EstProb[iii]
Dtct = Expect_Prob$Detect[iii]
SiteN = Expect_Prob$site[iii]

dfCallsL = data.frame(SiteN=as.integer(SiteN),WeekN=as.integer(WeekN),DayN = as.integer(DayN),
                      TS=as.integer(TS),Calls=Dtct,Moon = Moon)
dfCalls <- ddply(dfCallsL,  c("SiteN","WeekN","DayN","TS"), summarise,
               N    = length(Calls),
               Calls = sum(Calls),
               Moon=mean(Moon))
dfCalls$Minutes <- dfCalls$N/30 

dfProbL = data.frame(SiteN=as.integer(SiteN),WeekN=as.integer(WeekN),DayN = as.integer(DayN),
                     TS=as.integer(TS),probs=Prob,Moon = Moon, DNN = DNN, Flux=Flux,
                     LvAb = LvAb, Clck = Clck, Brst=Brst,Dtct = Dtct)
dfProb <- ddply(dfProbL,  c("SiteN","WeekN","DayN","TS"), summarise,
                Moon=mean(Moon),
                N    = length(probs),
                Prob = mean(probs),
                DNN = mean(DNN),
                Dtct = mean(Dtct),
                Flux = mean(Flux),
                LvAb = mean(LvAb),
                Clck = mean(Clck),
                Brst = mean(Brst))
dfProb$Flux2 <- dfProb$Flux^2
dfProb$LgtProb <- logit(dfProb$Prob)
dfProb$LgtDNN <- logit(pmax(.000001,pmin(.999999,dfProb$DNN)))
dfProb$Minutes <- dfProb$N/30 
dfProb$Calls <- dfProb$Dtct*dfProb$N

mod1 <- glm(Calls ~ DNN + Flux + LvAb + Clck + Brst + as.factor(SiteN), 
             data=dfProb, family=poisson(link = "log"))
summary(mod1)

mod2 <- lm(LgtDNN ~ Dtct + Flux + Flux2 + LvAb + Clck + Brst + as.factor(SiteN), 
            data=dfProb) #,family = gaussian(link="identity")
summary(mod2)


NObs = dim(dfCalls)[1]
NTsteps = max(dfCalls$TS)
TmprlMax = numeric(length=100)+6 # Ensures temporal matrix is 0-1

# Part 1: estimation using DNN and sub-set Audited ------------------------------------------------
AudPsite = 300
iii = numeric()
for (i in 1:NSite){
  iii = c(iii,sample(which(dfProb$SiteN==i),AudPsite,replace = FALSE))
}
NObsA = length(iii)
#
data <- list(LgtDNN=dfProb$LgtDNN,SiteN=dfProb$SiteN,N=dfProb$N,NSite=NSite,
             Minutes=dfProb$Minutes[iii],SiteNA=dfProb$SiteN[iii],Neff = 5,
             Flux=dfProb$Flux,Flux2=dfProb$Flux2,LvAb=dfProb$LvAb,
             Clck=dfProb$Clck,Brst=dfProb$Brst,
             Wk=dfProb$WeekN,TS=dfProb$TS,NObs=NObs,Moon=dfCalls$Moon, 
             WkA=dfProb$WeekN[iii],TSA=dfProb$TS[iii],Calls=dfProb$Calls[iii],
             NWeeks=NWeeks,NTsteps=NTsteps,NObsA=NObsA, 
             MoonA=dfCalls$Moon[iii], TmprlMax=TmprlMax) 
#
# Inits: Best to generate initial values using function
inits <- function(){
  list(sigT=runif(1,0.3,0.6),sigS=runif(1,2,10),sigA=runif(1,1,3),Dispers=runif(1,20,30))
}
# List of parameters to monitor:
params <- c('sigT','sigS','sigA','sigE','Dispers','phi','alph0','alpha','theta',
            'B','C','Cs','Temp','eps') # 
# Model to be run:
modfile = 'Jags_calls_mod3.jags'
# JAGS set-up parameters 
set.seed(123)
Nchains = 20
Nburnin =  5000  # Number of burn-in reps Total reps = (Nsim-Nburnin) * (num Cores)
Nadapt =  100  # Number of adapting reps, default 100
Totalreps = 10000 # Total desired reps (ie # simulations making up posterior)
Nsim =  Totalreps/Nchains + Nburnin 
# For parallel (comment out for serial)
cores<-detectCores()
cores = min(cores, Nchains)
cl <- makeCluster(cores[1])
registerDoParallel(cl)
#
out <- jags.basic(data = data,
                  inits = inits,
                  parameters.to.save = params,
                  model.file = modfile,
                  n.chains = Nchains,
                  n.adapt = Nadapt,
                  n.iter = Nsim,
                  n.burnin = Nburnin,
                  n.thin = 1,
                  parallel=TRUE,
                  n.cores=cores)

vn = varnames(out)
outmat = as.matrix(out)
reps = dim(outmat); reps = reps[1]
outdf = as.data.frame(outmat); rm(outmat)
s = summary(out) 
s_stats = data.frame(s$statistics)
s_quantiles = data.frame(s$quantiles) 

pfp = c(which(params=='theta'),which(startsWith(params,'sig')),which(params=='Dispers'),
        which(params=='C'),which(startsWith(params,'phi')),which(startsWith(params,'alph')))
for (i in pfp){
  parnm = params[i]
  traplot(out,parnm)
  denplot(out,parnm,ci=.9,collapse = TRUE)
}

pfp = c(which(startsWith(params,'Cs')))
for (i in pfp){
  parnm = params[i]
  traplot(out,parnm)
  denplot(out,parnm,ci=.9,collapse = TRUE)
}

# Plot heat map of the Temporal effect:
# Generate the "Temporal Effect" matrix from mean values
TemporalMn = matrix(0,nrow=NWeeks,ncol=NTsteps)
for (i in 1:NWeeks) {
  for (j in 1:NTsteps) {
    # TemporalMn[i,j] = mean(outdf[,which(vn==paste0('Temporal[',i,',',j,']'))])
    TemporalMn[i,j] = s_stats[which(vn==paste0('Temp[',i,',',j,']')),"Mean"]
  }
}
NDays = max(df$DayOfYear)-min(df$DayOfYear)
TS1 = parse_date_time('7:00PM',c('%I:%M %p'))
ii = which(df$DayOfYear==min(df$DayOfYear))
Date1 = ymd(df$Date[ii[1]])
Ntticks = 8
TS_ticks = signif(seq(1,NTsteps,len=Ntticks),3)
timelabels = as.character(format(TS1, format="%H:%M"))
for (t in 2:Ntticks){
  timelabels = c(timelabels,
                 as.character(format(TS1+(TS_ticks[t]-1)*15*60, format="%H:%M")))
}
Startdate = Date1
datevec = seq(Startdate+7, Startdate+NDays-1, by="2 weeks")
datelabels = format.Date(datevec,"%b")
# FortNtlabels = c('')

# Plot Image with color bar
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5,1), heights=c(1,1))
par(mar = c(3,5,2.5,2))
image(seq(1,NWeeks),seq(1,NTsteps),TemporalMn, col = heat.colors(100),
      xlab="Month",ylab="", yaxt="n", xaxt="n", cex.lab=.8)
title(main = "Temporal Effects on Call Rate", font.main = 4)
# add categorical labels to y-axis
axis(1, at=seq(2,NWeeks,by=2), labels=datevec, las=HORIZONTAL<-1, cex.axis=0.6)
axis(2, at=TS_ticks, labels=timelabels, las=HORIZONTAL<-1, cex.axis=0.6)
title(ylab = "Time of Day", line = 2.5, cex.lab=.8)

ColorLevels <- seq(0, 1, length=length(heat.colors(100)))
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col = heat.colors(100),
      xlab="",ylab="",
      xaxt="n", cex.axis=0.8, cex.lab=.8)
title(main = "Color Scale", font.main = 1)
title(xlab = "Proportion max rate", line = .5, cex.lab=.8)
layout(1)

# Part 2: estimation using Audited Calls ------------------------------------------------
data <- list(Calls=dfCalls$Calls,SiteN=dfCalls$SiteN,Minutes=dfCalls$Minutes,
             Wk=dfCalls$WeekN,TS=dfCalls$TS,NObs=NObs,NSite=NSite, 
             NWeeks=NWeeks,NTsteps=NTsteps,
             Moon=dfCalls$Moon,TmprlMax=TmprlMax) 

# Inits: Best to generate initial values using function
inits <- function(){
  list(sigT=runif(1,0.3,0.6),sigS=runif(1,2,10),sigW=runif(1,.1,.5),
       Dispers=runif(1,.2,.4))
}
# List of parameters to monitor:
paramsA <- c('theta','sigT','sigS','sigW','Dispers',
            'C','Cs','Temp','eps') # 
#     NOTE: to save "E" matrix, need to run in serial not parallel
#           due to memory limitations of parallel

# Model to be run:
modfile = 'Jags_calls_mod1.jags'
# JAGS set-up parameters 
set.seed(123)
Nchains = 20
Nburnin =  5000  # Number of burn-in reps Total reps = (Nsim-Nburnin) * (num Cores)
Nadapt =  100  # Number of adapting reps, default 100
Totalreps = 10000 # Total desired reps (ie # simulations making up posterior)
Nsim =  Totalreps/Nchains + Nburnin 
# For parallel (comment out for serial)
cores<-detectCores()
cores = min(cores, Nchains)
cl <- makeCluster(cores[1])
registerDoParallel(cl)
#
outA <- jags.basic(data = data,
                  inits = inits,
                  parameters.to.save = paramsA,
                  model.file = modfile,
                  n.chains = Nchains,
                  n.adapt = Nadapt,
                  n.iter = Nsim,
                  n.burnin = Nburnin,
                  n.thin = 1,
                  parallel=TRUE,
                  n.cores=cores)

vnA = varnames(outA)
outmatA = as.matrix(outA)
repsA = dim(outmatA); repsA = repsA[1]
outdfA = as.data.frame(outmatA); rm(outmatA)
sA = summary(outA) 
sA_stats = data.frame(sA$statistics)
sA_quantiles = data.frame(sA$quantiles) 

pfp = c(which(paramsA=='theta'),which(startsWith(paramsA,'sig')),which(paramsA=='Dispers'),
        which(paramsA=='C'),which(paramsA=='psi'))
for (i in pfp){
  parnm = paramsA[i]
  traplot(outA,parnm)
  denplot(outA,parnm,ci=.9,collapse = TRUE)
}
pfp = c(which(startsWith(paramsA,'Cs')))
for (i in pfp){
  parnm = paramsA[i]
  traplot(outA,parnm)
  denplot(outA,parnm,ci=.9,collapse = TRUE)
}

# Part 3: estimation using Estiomated Probs only ------------------------------------------------
data <- list(LgtProb=dfProb$LgtProb,SiteN=dfProb$SiteN,N=dfProb$N,
             Wk=dfProb$WeekN,TS=dfProb$TS,NObs=NObs,NSite=NSite, 
             NWeeks=NWeeks,NTsteps=NTsteps,
             Moon=dfCalls$Moon,TmprlMax=TmprlMax) 

# Inits: Best to generate initial values using function
inits <- function(){
  list(sigT=runif(1,0.3,0.6),sigS=runif(1,2,10),sigW=runif(1,.1,.5),
       sigE=runif(1,1,2))
}
# List of parameters to monitor:
paramsP <- c('theta','sigT','sigS','sigW','sigE',
            'C','Cs','Temp','eps') # 
#     NOTE: to save "E" matrix, need to run in serial not parallel
#           due to memory limitations of parallel

# Model to be run:
modfile = 'Jags_calls_mod2.jags'
# JAGS set-up parameters 
set.seed(123)
Nchains = 20
Nburnin =  5000  # Number of burn-in reps Total reps = (Nsim-Nburnin) * (num Cores)
Nadapt =  100  # Number of adapting reps, default 100
Totalreps = 10000 # Total desired reps (ie # simulations making up posterior)
Nsim =  Totalreps/Nchains + Nburnin 
# For parallel (comment out for serial)
cores<-detectCores()
cores = min(cores, Nchains)
# cl <- makeCluster(cores[1])
registerDoParallel(cl)
#
outP <- jags.basic(data = data,
                  inits = inits,
                  parameters.to.save = paramsP,
                  model.file = modfile,
                  n.chains = Nchains,
                  n.adapt = Nadapt,
                  n.iter = Nsim,
                  n.burnin = Nburnin,
                  n.thin = 1,
                  parallel=TRUE,
                  n.cores=cores)

vnP = varnames(outP)
outmatP = as.matrix(outP)
repsP = dim(outmatP); repsP = repsP[1]
outdfP = as.data.frame(outmatP); rm(outmatP)
sP = summary(outP) 
sP_stats = data.frame(sP$statistics)
sP_quantiles = data.frame(sP$quantiles) 

pfp = c(which(paramsP=='theta'),which(startsWith(paramsP,'sig')),which(paramsP=='Dispers'),
        which(paramsP=='C'),which(paramsP=='psi'))
for (i in pfp){
  parnm = paramsP[i]
  traplot(outP,parnm)
  denplot(outP,parnm,ci=.9,collapse = TRUE)
}
pfp = c(which(startsWith(paramsP,'Cs')))
for (i in pfp){
  parnm = paramsP[i]
  traplot(outP,parnm)
  denplot(outP,parnm,ci=.9,collapse = TRUE)
}

save.image("./files/compare_Bayescalls_Audit_Prob.rdata")

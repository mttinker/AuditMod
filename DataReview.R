# Data Review
df = read.csv("OurData_Subset_Audited.csv")
TargetCall = c("HAPE_h")
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
tmp = unique(df$rating[Detect==0 & Prob1>.95])

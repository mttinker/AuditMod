# Shell to run Function to fit probability model to audited samples for varying sample sizes 
df = read.csv("OurData_Subset_Audited.csv")
TargetCall = c("HAPE_h")
source("DataReview.r")
rslt = DataReview(df,TargetCall)
load(rslt)
rm(DataReview)
AuditSzs = c(100,300,500,750,1000)
NSS = length(AuditSzs)
# Text1 = c("./files/fitAuditProb_Result_7s_")
# Text1b = c("./files/fitAuditProb_Small_7s_")
# Text2 = c("_AudPerSite.rdata")
# FileList = character(length = NSS)
# Savenames = character(length = NSS)
# for (f in 1:NSS){
#   Savenames[f] = paste0(Text1,AuditSzs[f],Text2)
#   FileList[f] = paste0(Text1b,AuditSzs[f],Text2)
# }
source("FitAuditProb.r")
rslts = vector("list", NSS)
# Run Bayesian model to fit Probability to audited data for varying sample sizes---------------
#
for (k in 1:NSS){
  rslt = FitAuditProb(AuditSzs[k],df,Sitelist)
  rslts[k] = rslt
}
source("Audit_Sample_Power.r")
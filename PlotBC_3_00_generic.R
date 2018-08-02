# Set-up for doing plots from BayesCalls analysis
# Purpose of this file is for user to specify parameters needed for plotting results
# NOTE: Assumed that this script is in same directory as Results file
rm(list = ls())
# USER SPECIFIED PARAMETERS -------------------------------------------------
Resultsfile = "Results_HAPE_2017_BC3_00_2018-08-02-13h"
ResultsFolder = 'AuditMod/files'  # Folder path within RootDir where results files stored
# Set root directory path... enter absolute path or relative,
RootDir =  "above1"  # Examples "current" or "above1" or "C:/Users/XXXX/BayesianMethods"
AnalysisFolder = 'Acoustic3'  # Folder path within RootDir where plot code is stored
DataFolder = 'AuditMod/files'  # Folder path within RootDir where raw data files stored
RunFile = 'BayesCallsPlots3_00'       # Plotting Script
# Next lines used to select sites for making plots of mean expected call rate
#  over a focal period (can be different from peak period)
Nfocalsites = 3 # Number of sites defined below... generally 10 or less
Site_focal = c("NBOG4","NBOG5","ULP8")
# Define focal period of interest
DateFocalStrt = "June 1"
DateFocalStop = "June 30"
TimeFocalStrt = 60  # Time boundary 1, minutes relative to a reference event (sunrise or sunset)
TimeFocalStop = 120  # Time boundary 2, always > than boundary1 (minutes relative to event) 
Time_ref = 1  # Reference event: 1 = after sunset, 2 = before sunrise
# Specify name of Data files with Conversion function and Areas of each strata
Convertfile = ''
Areasdatfile = '' # Name of matching data file with nest count data (OTIONAL, enter blank if no nest counts)

# END USER SPECIFIED PARAMETERS ---------------------------------------------
#
# (Should not need to edit anything below here)
# 
# Process Filenames and directory names -------------------------------------
#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
thisdir = getwd()
if (RootDir=='current'){
  RootDir = getwd()
} else if(RootDir=='above1'){
  tmp1 = getwd()
  setwd('..')
  RootDir = getwd()
  setwd(tmp1)
  rm('tmp1')
} else if(RootDir=='above2'){
  tmp1 = getwd()
  setwd('..')
  setwd('..')
  RootDir = getwd()
  setwd(tmp1)
  rm('tmp1')
} else {
  RootDir = RootDir
}

RunFile = paste0(RootDir,"/",AnalysisFolder,"/",RunFile,".r")
ResultsFolder = paste0(RootDir,"/",ResultsFolder)
DataFolder = paste0(RootDir,"/",DataFolder)
Areasdatfile = paste0(DataFolder,"/",Areasdatfile)
loadfile1 = paste0(ResultsFolder,"/",Resultsfile,".rdata")
loadfile2 = paste0(ResultsFolder,"/",Convertfile,".rdata")
AnalysisFolder = paste0(RootDir,"/",AnalysisFolder)
rm(thisdir)
#
source(RunFile)

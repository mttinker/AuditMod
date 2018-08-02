# RunBC_EDIT_THIS_PART = User File to Run "BayesCalls" Program to analyze Acoustic call data 
#    using a Bayesian State Space model (BayesCalls)
#    (***filename needs to begin with "RunBC_")
# Purpose of this file is for user to specify parameters needed to run analysis
# NOTE: Results file will get same name as this script, with "RunBC" replaced by "Results"
rm(list = ls())
gc()
existing_Packages<-as.list(installed.packages()[,1])
# Add new packages you might need to this line, to check if they're installed and install missing packages
required_Packages<-c("rstudioapi")
missing_Packages<- required_Packages[!required_Packages%in% existing_Packages]
if(length(missing_Packages)>0)install.packages(pkgs =  missing_Packages)
invisible(lapply(required_Packages, require, character.only=T,quietly = T))
rm(existing_Packages)
start.time <- Sys.time()  # start a timer
# USER SPECIFIED PARAMETERS -------------------------------------------------
#
# First provide a "Root Directory" within which all other folders, subfolders & files are contained
#   NOTE: enter "current" for current folder, "above1" for 1 folder above, "above2" for 2 above, 
#   otherwise enter path with "/" separator for sub-folders
RootDir =  "above1"  # Examples "current" or "above1" or "C:/Users/XXXX/BayesianMethods"
AnalysisFolder = 'Acoustic3'  # Folder path within RootDir where analysis code is stored
RunFile = 'BayesCalls3_00'       # Version of BayesCalls to run
DataFolder = 'AuditMod/files'  # Folder path within RootDir where raw data files stored
ResultsFolder = 'AuditMod/files'  # Folder path within RootDir where results files stored
ProjectName =  'Bayesian_Auditmod_2017'  # Name of project (for user reference)
Species =  'HAPE'  # Species name for data analysis
Countsdatfile = c("") # Name of matching data file with nest count data (OTIONAL, enter blank if no nest counts)
Datafile =  paste0('HAPE_15min_all_prob_1k1k.RData')  #Name of data file for analysis 
ProjectLocation =  'Kuai'  # Island name for data analysis
ProjectYear =  2017  # Year of analysis project 
Yearfocal =  2017  # Focal year for Bayesian analysis
TS_strt =  1  # Starting time step (beginning of range to load data)
TS_stop =  46  # Ending timestep (end of range to load data)
TimeStep1 = '18:45'  # Clock time (24hr) of timestep1
TimeStepIntvl =  15  # Number of minutes in each timestep
subsamp =  3  # Level of Sub-sampling of entire data set: use every nth record
DayOfYear_strt =  100  # April 10
DayOfYear_end =  206  # July 25
peakdates_strt = 125 # May 5
peakdates_end = 155  # June 4
peaktimes_strt =  60  # Peak time boundary 1, minutes relative to a reference event (sunrise or sunset)
peaktimes_stop =  120  # Peak time boundary 2, always > than boundary1 (minutes relative to event) 
peaktimes_ref =  1  # Reference event: 1 = after sunset, 2 = before sunrise, 3 = sunrise AND sunset 
calendar_opt = 1  # year range: 1 = All data within calendar year, 2 = data spans New Year
calendar_pk_opt = 1 # 1 = all in one year, 2 = peak spans new year
data_opt =  1  # Data Option: 1 = Calls Only, 2 = Calls plus Nest Counts
prior_opt = 1   # Priors: 1 = uninformed, 2 = informed (must supply results file)
PriorResultsFile =  c()  # OPTIONAL: if prior_opt = 2, Rdata file containing parameter priors
#
Nchains = 20
Nburnin =  5000  # Number of burn-in reps Total reps = (Nsim-Nburnin) * (num Cores)
Nadapt =  500  # Number of adapting reps, default 100
Totalreps = 5000 # Total desired reps (ie # simulations making up posterior)
#
# END USER SPECIFIED PARAMETERS ---------------------------------------------
#
# (Should not need to edit anything below here)
# 
Nsim =  Totalreps/Nchains + Nburnin  # Total # MCMS sims: Actual saved reps = (Nsim-Nburnin) * (num Cores)
# Process Filenames and directory names -------------------------------------
#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
RunFile1 = RunFile
RunFile1 = gsub("BayesCalls","BC",RunFile1)
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
#
if (AnalysisFolder=='') {
    RunFile = paste0(RootDir,"/",RunFile,".R")
    AnalysisFolder = getwd()
} else {
    RunFile = paste0(RootDir,"/",AnalysisFolder,"/",RunFile,".R")
    AnalysisFolder = paste0(RootDir,"/",AnalysisFolder)
}
#
if (DataFolder=='') {
    loadfile1 = paste0(RootDir,"/",Datafile)
} else {
    loadfile1 = paste0(RootDir,"/",DataFolder,"/",Datafile)
}
#
if (data_opt>1){
  if (DataFolder=='') {
    loadfile2 =  paste0(RootDir,"/",Countsdatfile)
  } else {
    loadfile2 =  paste0(RootDir,"/",DataFolder,"/",Countsdatfile)
  }
}

DataFolder = paste0(RootDir,"/",DataFolder)
ResultsFolder = paste0(RootDir,"/",ResultsFolder)
#
#  Function to get the name of this current user script 
callscript <- function() {
    # http://stackoverflow.com/a/32016824/2292993
    cmdArgs = commandArgs(trailingOnly = FALSE)
    needle = "--file="
    match = grep(needle, cmdArgs)
    if (length(match) > 0) {
        # Rscript via command line
        return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
        ls_vars = ls(sys.frames()[[1]])
        if ("fileName" %in% ls_vars) {
            # Source'd via RStudio
            return(normalizePath(sys.frames()[[1]]$fileName)) 
        } else {
            if (!is.null(sys.frames()[[1]]$ofile)) {
                # Source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
            } else {
                # RStudio Run Selection
                # http://stackoverflow.com/a/35842176/2292993  
                return(normalizePath(rstudioapi::getActiveDocumentContext()$path))
            }
        }
    }
}
# Get the name of current script file, and create matching SaveResults filename 
CallingScript = callscript()
SaveResults = gsub("RunBC","Results",CallingScript)
SaveResults = paste0(SaveResults,"data")
# Alternative: Standard Save results name including species, year, runfile name and date
SaveResults1 = RunFile1
tmp = as.character(start.time)
SaveResults1 = paste0(ResultsFolder,'/Results', '_', Species, '_', ProjectYear, '_', 
                      SaveResults1,"_",substr(tmp,1,10),"-",substr(tmp,12,13),"h.Rdata")
SaveResults = SaveResults1
rm('tmp')
rm('SaveResults1')
rm('callscript')
#
# Run Bayesian Analysis -----------------------------------------------------
# Source BayesCalls program (specify version as X_XX, e.g. BayesCalls2_01) 
#
source(RunFile)

Totalend.time <- Sys.time()
Totaltime.taken <- start.time - Totalend.time
Totaltime.taken


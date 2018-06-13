# CMI GENERIC R Code to extract data from CMI_Auditor log files

# Intent: This script takes log files from CMI_Auditor and creates our
# analysis_data file which is saved as an RDS file.  This includes taking 2
# second window s and wrapping them up by minute to have a count per minute of
# the signal. We also take "p" rows which there will be one row per file (if
# files are <1 minute) or 1 row per minute (if files are >1 minute) and fills in
# zeros if there are not detections during that minute.  This allows us to
# quantify recording effort. We also add field metadata to the data as well as
# sunrise/set info
#
# Input: log files, field data csv
# Outputs: Analysis data .Rda file
# 
# Limitations: currently the Almanac function needs to have internet connectivity 
# and the navy website needs to be up and running


rm(list=ls())
start.time <- Sys.time()  # start a timer

# Input variables -----------------------------------------------------------------------------

pathToRepo="D:/CM,Inc/R-Code_improvements"
Dropbox="D:/CM,Inc/Dropbox (CMI)"

ProjectName="Test_Dataset_2018"                             # Name of project folder from Dropbox
ProjectYear=2018                                            # Year of project from Dropbox
Species="HAPE"                                              # Species code to name outputs with
ProjectLocation="Kauai"                                     # Island name
deployment_info<-""              # name of field data csv
Round=""                                                 # round name to name outputs with
Nocturnal=T                                                 # True = nocturnal species (attach sunrise/set info by night)
RemoveIntro=F                                               # do you want to remove intro recordings?

proj4string <- "+proj=utm +zone=4 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # must change to correct zone and north or south of equator
Almanac_Input<-list(Years=c("2017"), # Year or vector of years
                    LatD="22", # Degrees Latitude of center of project area (MUST be two digits and quoted)
                    LatM="09", # Minutes Latitude of center of project area (MUST be two digits and quoted)
                    LonD="159", # Degrees Logitude of center of project area (MUST be two digits and quoted)
                    LonM="36", # Minutes Longitude of center of project area (MUST be two digits and quoted)
                    East_West="-1", # -1 for minus Longitude (west of prime maridian)
                    North_South="1", #-1 for south of equator
                    TimeZoneName = "US/Hawaii", # this is the name of the time zone find at use OlsonNames() to get a list
                    TimeZ="10", #time zone from GMT (MUST be two digits and quoted)
                    TimeZ_East_West="-1", # -1 for minus GMT (west of GMT)
                    Location=ProjectLocation, # Do Not Change  
                    Hour="00" # for moon illumination: the hour of the day.  00 is midnight
)

if (Species=='HAPE') {
    saveTagSubsetAudit="_Audit_HAPE_HAPE_21Apr"
    saveTagFull="audit_HAPE_21Apr"
    LabelPositive='HAPE_h' # positive label in the fig, default Species code
    threshold=0.95
} else if (Species=='NESH') {
    saveTag='audit_NESH_16Apr'
}

# Prep workspace ----------------------------------------------------------

ProjectDir='D:/CM,Inc/Dropbox (CMI)/CMI_Team/Analysis/Audit_Subset_Workflow'

source('D:/CM,Inc/R-Code_improvements/Functions/Almanac_Maker2.R')

existing_Packages<-as.list(installed.packages()[,1])
required_Packages<-c("data.table", "dplyr", "tidyr", "stringr", "readr", "lubridate",
                     "httr","zoo", "oce", "ggplot2","gridExtra",
                     "proj4", "sf", "mapview","rgdal",
                     "openxlsx")
missing_Packages<- required_Packages[!required_Packages  %in% existing_Packages]
if(length(missing_Packages)>0)install.packages(pkgs =  missing_Packages)
print(paste("Loading Packages: ",paste(required_Packages,collapse =", " )))
invisible(suppressMessages(lapply(required_Packages, require, character.only=T,quietly = T)))

# Load Data -------------------------------------------------------------------

dirs<-list.dirs(paste(ProjectDir,"Log_Files",Round,sep="/"))
# audited subset figs
file.list.subset <- list.files(dirs, pattern=paste0(saveTagSubsetAudit,".csv"), full.names=T, recursive=F)
if (length(file.list.subset)>15) {head(file.list.subset,10)} else {file.list.subset}
# fully audited figs > threshold prob
file.list.full <- list.files(dirs, pattern=paste0(saveTagFull,"(_KJD)?.csv"), full.names=T, recursive=F)
if (length(file.list.full)>15) {head(file.list.full,10)} else {file.list.full}

# Read log files for subset audited figs
data.list.subset <- lapply(file.list.subset, fread, sep=",",skip=" ",showProgress=T, drop=c( "notes","date", "time"))
head(data.list.subset,1)
length(file.list.subset) # double check that this is how many figs you should have
# REad log files for fully audited figs above threshold
data.list.full <- lapply(file.list.full, fread, sep=",",skip=" ",showProgress=T, drop=c("notes","date", "time"))
head(data.list.full,1)
length(file.list.full) # double check that this is how many figs you should have

# Make OurData -------------------------------------------

# # make audit_file column if missing
# for(i in c(1:length(data.list.subset))){
#     if(!'audit_file' %in% names(data.list.subset[[i]])){
#         data.list.subset[[i]]$audit_file<-NA
#     }
# }

# concatenate into two big data.frames
OurDataSubset <- as.data.frame(bind_rows(data.list.subset)) %>% filter(rating!='P')
head(OurDataSubset)
OurDataFull <- as.data.frame(bind_rows(data.list.full)) %>% filter(rating!='P')
head(OurDataFull)

# filter out weird repeated rows
OurDataSubset <- OurDataSubset %>% filter(!is.na(id))

# check for NA ratings
table(OurDataSubset$rating)
sum(is.na(OurDataSubset$rating))
table(OurDataFull$rating)
sum(is.na(OurDataFull$rating))

sapply(OurDataSubset,function(x) sum(is.na(x)))
sapply(OurDataFull,function(x) sum(is.na(x)))

OurDataSubset$key<-paste(OurDataSubset$file,OurDataSubset$start_in_file)
OurDataFull$key<-paste(OurDataFull$file,OurDataFull$start_in_file)

length(which(!(OurDataSubset$key %in% OurDataFull$key)))
table(OurDataSubset$location[!(OurDataSubset$key %in% OurDataFull$key)])

OurData<-OurDataSubset %>% rename('rating_new'='rating') %>% left_join(select(OurDataFull,key,'rating_old'='rating'),by=c('key'))

table(OurData$rating_new)
table(OurData$rating_new,OurData$rating_old)
sum(is.na(OurData$rating_new))
sum(is.na(OurData$rating_old)) # ok to have NAs here bc we only auditing >threshold before
table(OurData$rating_new[is.na(OurData$rating_old)]) # make sure we have new ratings for old NAs
table(OurData$location[is.na(OurData$rating_old) & OurData$rating_new=='HAPE_hi']) # make sure we have new ratings for old NAs

OurData$rating=ifelse(OurData$rating_new=='HAPE_hi',OurData$rating_old,OurData$rating_new)
table(OurData$rating)
sum(is.na(OurData$rating))
table(OurData$rating_new[is.na(OurData$rating)]) # make sure we have new ratings for old NAs
table(OurData$rating_old[is.na(OurData$rating)]) # make sure we have new ratings for old NAs
sum(is.na(OurData$rating_old[is.na(OurData$rating)])) # make sure we have new ratings for old NAs
OurData$rating[OurData$rating=='4']='HAPE_tone'
OurData$rating[OurData$rating=='5']='HAPE_h'
table(OurData$rating)

unique(OurData$audit_file[OurData$rating=='HAPE_mid'])

# Fix datetimes ---------------------------------------------------------------

# make initial DateTime
OurData$DateTime<-ymd_hms(paste(OurData$year, OurData$month, OurData$day,
                                OurData$hour,OurData$minute,OurData$second),
                          tz=Almanac_Input$TimeZoneName)

# Find start second of recordings
OurData$rec_start_second<-str_extract(OurData$file,"(?<!\\d)\\d{6}(?!\\d)")
OurData$rec_start_second<-as.numeric(str_extract(OurData$rec_start_second,"[0-9]{2}$"))
head(OurData)

# subtract off the start seconds so that all the date times are shifted as if
# the recordings started at 00 seconds
OurData$DateTime<-OurData$DateTime-OurData$rec_start_second

# repopulate all date and time objects
OurData$year<-year(OurData$DateTime)
OurData$month<-month(OurData$DateTime)
OurData$day<-day(OurData$DateTime)
OurData$hour<-hour(OurData$DateTime)
OurData$minute<-minute(OurData$DateTime)
OurData$second<-second(OurData$DateTime)
OurData$Date<-as.Date(OurData$DateTime)

# Remove Metrics ("P" values) to be added back in later -----------------------
# 
# # Make probeId key
# OurData$key<-paste(OurData$audit_file,OurData$probeId,sep="_")
# 
# if(("P" %in% OurData$rating)){
#     print("There are entries with rating==P")
#     Metrics<-subset(OurData,rating=="P")
#     OurData<-subset(OurData,rating!="P")
# }else{print("There are NO entries with rating==P")}
# 
# # check for duplicated probeId/audit_file keys
# ifelse(sum(duplicated(Metrics$key))==0,"You're G as can be!","You have duplicated probeId keys--You probably pulled in duplicate log files.")

# Count Hits ----------------------------------------------------------------------------------

# # Rows with Rating=5 (or species 4-letter code) get a 1 while all others get 0
# OurData$Total_hits<-ifelse(OurData$rating==LabelPositive,1,0)
# OurData$Total_alarms<-ifelse(OurData$rating!=LabelPositive,1,0)
# # Fill in any NAs
# OurData$Total_hits[is.na(OurData$Total_hits)]<-0
# OurData$Total_alarms[is.na(OurData$Total_alarms)]<-0

# Data Check 1 --------------------------------------------------------------------------------

# table(OurData$rating)  
# unique(OurData$Total_hits)
# unique(OurData$Total_alarms)
# 
# DataCheck1<-OurData %>%
#     group_by(location) %>%
#     summarise(OurDataTotals=sum(Total_hits))
# as.data.frame(DataCheck1)

# summarise by minute -------------------------------------------------------------------------

# Build Date
OurData$Date<-ymd(paste(OurData$year, OurData$month, OurData$day))

# # Use dplyr to summarise
# HitsperMin<-OurData %>%
#     group_by(key) %>%
#     summarise(Total_hits=sum(Total_hits),Total_alarms=sum(Total_alarms))

# Add Metrics to data -------------------------------------------------------------------------

# # This also adds any minutes of effort with no detections back into the data with 0s
# ifelse(sum(duplicated(Metrics$key))==0,"Still no duplicated P rows! Continue.","Uh oh! You have duplicated P rows!")
# Metrics$Date<-ymd(paste(Metrics$year, Metrics$month, Metrics$day))
# 
# OurData1<-left_join(Metrics,HitsperMin,by="key")
# ifelse(nrow(Metrics)==nrow(OurData1),"Everything joined correctly!","Something went wrong with the join...Maybe you had duplicated keys?")
# 
# # Populate NAs for minutes that had no detections ---------------------------------------------
# 
# OurData1$Total_hits[is.na(OurData1$Total_hits)]<-0
# OurData1$Total_alarms[is.na(OurData1$Total_alarms)]<-0

# Add date time and jdate ---------------------------------------------------------------------

OurData$DateTime<-ymd_hm(paste(OurData$year, OurData$month, OurData$day, OurData$hour, OurData$minute),tz=Almanac_Input$TimeZoneName)
OurData$DayOfYear<-as.numeric(lubridate::yday(OurData$DateTime))

# # Data Check 2 --------------------------------------------------------------------------------
# 
# DataCheck2<-OurData%>%
#     group_by(location)%>%
#     summarise(OurDataTotals=sum(Total_hits))
# DataCheck<-full_join(DataCheck1,DataCheck2)
# DataCheck$diff<-DataCheck$OurDataTotals-DataCheck$OurDataTotals
# if(sum(DataCheck$diff)!=0) stop("You are loosing or gaining strikes")
# as.data.frame(DataCheck)

# Add SPID ------------------------------------------------------------------------------------

# # read in SPID field data
# deployPath<-list.files(ProjectDir,pattern=deployment_info,full.names=T, recursive=T)
# SPdata <- read.csv(deployPath,stringsAsFactors = F)
# # Work around if search doesn't work
# # SPdata <- read.csv(file.path(ProjectDir,'Field_Data','SM_Deployment_Info_KESRP_ULP_HON_2013_all_rounds_FixedSensor_Name_AF.csv'))
# 
# SPdata<-spDataCleanR(SPdata)

# Add SPID  --------------------------------------------------------------------

# # Make it a data.table
# setDT(SPdata)
# 
# #set the index for the SPdata whick is what allows to magic in foverlaps
# setkey(SPdata, Sensor_Name,Deployment_DateTime, Retrieval_DateTime)
# 
# OurData$Sensor_Name<-OurData$location
# OurData$DateTime2<-OurData$DateTime
# 
# # Make it a data.table
# setDT(OurData)
# 
# # fast function from data.tble that joins by a time interval
# OurData<-foverlaps(x = OurData,y = SPdata,
#                          by.x=c("Sensor_Name","DateTime","DateTime2"),
#                          mult="all",type="any") %>% 
#     filter(!is.na(SPID)) %>%
#     arrange(SPID,DateTime)
# 
# sum(duplicated(OurData$key)) # should be 0

# Make missing effort graph ---------------------------------------------------------------------------------

# # data that gets removed by SPdata
# NoSPID<-filter(OurData,!key %in% OurData$key)
# nrow(NoSPID) # how many minutes get cut out?
# 
# if (nrow(NoSPID) > 0) {
#     
#     unique(NoSPID$Sensor_Name)
#     
#     OurData$type=ifelse(OurData$key %in% OurData$key,'corrected','lost')
#     MissingEffortFig=makeEffort(AD=OurData,siteName='Sensor_Name',colors=c('red','grey70'))
#     
#     # save missing effort figure(s)
#     for (i in 1:(length(MissingEffortFig)-1)) {
#         pathMissingEffortFig<-file.path(ProjectDir,"Field_Data",paste0("Effort_Missing_",i,'_',Today,".jpg"))
#         ggsave(MissingEffortFig[[i]], filename=pathMissingEffortFig, width=6.5, height=8.4, units='in', scale=1)
#     }
#     
#     # save summary csv of missing effort
#     write.csv(MissingEffortFig[[length(MissingEffortFig)]],file.path(ProjectDir,"Field_Data",paste0("Effort_Missing_",Today,".csv")),row.names=F)
# 
# }

# Add NightOfYear -----------------------------------------------------------------------------

OurData$NightOfYear<-ifelse(hour(OurData$DateTime)%in%12:23,yday(OurData$Date),(yday(OurData$Date)-1))
OurData$SRiseNight<-ifelse(hour(OurData$DateTime)%in%0:11,yday(OurData$Date),(yday(OurData$Date)+1))
OurData$NightDate<-as.Date(OurData$NightOfYear, origin=paste(OurData$year,"-01-01",sep=""))- days(1)
OurData$SRiseNightDate<-as.Date(OurData$SRiseNight, origin=paste(OurData$year,"-01-01",sep=""))- days(1)

# if you pass from one year to the next you must use these lines so that day of year does not wrap. 
if(length(unique(OurData$year))==2){
    OurData$NightOfSeason<-ifelse(year(OurData$DateTime)==min(OurData$year),
        OurData$NightOfYear,OurData$NightOfYear+getDays(year = min(OurData$year)))
}
head(OurData)

# Data Check 3 --------------------------------------------------------------------------------

# DataCheck3<-OurData%>%
#     group_by(location)%>%
#     summarise(OurDataTotals=sum(Total_hits))
# 
# DataCheck<-full_join(DataCheck,DataCheck3)
# DataCheck$diff2<-DataCheck$OurDataTotals-DataCheck$OurDataTotals
# if(sum(DataCheck$diff2)!=0) stop("You are loosing or gaining detections")
# as.data.frame(DataCheck)

# Get Almanac ---------------------------------------------------------------------------------

# make list of dates
AlmanacDates<-data.frame(Date=seq(as.Date(min(OurData$Date))-1,as.Date(max(OurData$Date))+1,by='1 day'))
AlmanacDates$DateTime<-as.POSIXct(paste0(AlmanacDates$Date,' ',Almanac_Input$Hour,':00:00'),tz=Almanac_Input$TimeZoneName)
# new Almanac_Maker2 function
Almanac<-Almanac_Maker2(Lon = ((as.numeric(Almanac_Input$LonD)+(as.numeric(Almanac_Input$LonM)/60))*as.numeric(Almanac_Input$East_West)),
                        Lat = ((as.numeric(Almanac_Input$LatD)+(as.numeric(Almanac_Input$LatM)/60))*as.numeric(Almanac_Input$North_South)),
                        DateTime=AlmanacDates$DateTime,
                        tz=Almanac_Input$TimeZoneName)

# OR USE OLD ALMANAC
# Almanac<-Almanac_Maker(Almanac_Input$Years,
#                        Almanac_Input$LatD,
#                        Almanac_Input$LatM,
#                        Almanac_Input$LonD,
#                        Almanac_Input$LonM,
#                        Almanac_Input$East_West,
#                        Almanac_Input$North_South,
#                        Almanac_Input$TimeZ,
#                        Almanac_Input$TimeZ_East_West,
#                        Almanac_Input$Location,
#                        Almanac_Input$type,
#                        Almanac_Input$Hour)

# # OR BRING IN ALREADY MADE ALMANAC
# Almanac<-readRDS(file.path(Dropbox,"Analysis","Almanac_2012-2017_Kauai.Rda"))
# head(Almanac)

# Remove a few columns before you merge so that you do not have duplicate names
Almanac<-Almanac[,!(names(Almanac) %in%  c("Latitude","Longitude","Location"))]

if ('moonrise_time' %in% names(Almanac)) {
    AlmanacRise<-select(Almanac,Date,sunrise_time,moonrise_time,civil_dawn,naut_dawn,astronomical_dawn,Illu)
    AlmanacSet<-select(Almanac,-c(sunrise_time,moonrise_time,civil_dawn,naut_dawn,astronomical_dawn,Illu)) 
} else {
    AlmanacRise<-select(Almanac,Date,sunrise_time,civil_dawn,naut_dawn,astronomical_dawn,Illu)
    AlmanacSet<-select(Almanac,-c(sunrise_time,civil_dawn,naut_dawn,astronomical_dawn,Illu))
}

# Merge almanac with Combined data
if(Nocturnal) OurData<-left_join(OurData,AlmanacRise,by=c("SRiseNightDate" = "Date"))
if(!Nocturnal) OurData<-left_join(OurData,AlmanacRise,by=c("Date" = "Date"))

# Merge almanac with Combined data
if(Nocturnal) OurData<-left_join(OurData,AlmanacSet,by=c("NightDate" = "Date"))
if(!Nocturnal) OurData<-left_join(OurData,AlmanacSet,by=c("Date" = "Date"))

# Remove intro recordings ---------------------------------------------------------------------

# Remove the 3-min recordings from deployment (Sometimes 5 minutes)  
# Currently does not get ride of intro recordings that were not complete when the unit was turned off before the recording was finished
if(RemoveIntro){OurData<-subset(OurData,round(file_duration_sec)!=180|round(file_duration_sec)!=300) }

# Calculate DiffTimes -------------------------------------------------------------------------

# Calculate diff times for all sun/moon
OurData$MinFromSunset<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$sunset_time)))/60
OurData$MinFromSunrise<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$sunrise_time)))/60
if ('moonset_time' %in% names(OurData)) {
    OurData$MinFromMoonset<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$moonset_time)))/60
    OurData$MinFromMoonrise<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$moonrise_time)))/60
}
OurData$MinFromCivilDawn<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$civil_dawn)))/60
OurData$MinFromCivilDusk<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$civil_dusk)))/60
OurData$MinFromNautDawn<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$naut_dawn)))/60
OurData$MinFromNautDusk<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$naut_dusk)))/60
OurData$MinFromAstroDawn<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$astronomical_dawn)))/60
OurData$MinFromAstroDusk<-as.numeric((as.numeric(OurData$DateTime)-as.numeric(OurData$astronomical_dusk)))/60

#  Create your final AnalysisDataFull -------------------------------------------------------------------------------------------

OurData<-select(OurData,-rating_new,-rating_old)
names(OurData)

# # the below syntax should be robust to columns that are missing add anything that was on SPID table that you want in the final data
# names(OurData)
# AnalysisDataFull<-data.frame(SPID=OurData$SPID)
# 
# if(exists("SPID", OurData)) AnalysisDataFull$SPID<-OurData$SPID
# if(exists("location", OurData)) AnalysisDataFull$Sensor_Name<-OurData$location
# if(exists("Site_Name", OurData)) AnalysisDataFull$Site_Name<-OurData$Site_Name
# if(exists("Island", OurData)) AnalysisDataFull$Island<-OurData$Island
# if(exists("Region", OurData)) AnalysisDataFull$Region<-OurData$Region
# if(exists("SuperRegion", OurData)) AnalysisDataFull$SuperRegion<-OurData$SuperRegion
# if(exists("Site_Description", OurData)) AnalysisDataFull$Site_Description<-OurData$Site_Description
# if(exists("Site_Type", OurData)) AnalysisDataFull$Site_Type<-OurData$Site_Type
# if(exists("Predators", OurData)) AnalysisDataFull$Predators<-OurData$Predators
# if(exists("Elevation", OurData)) AnalysisDataFull$Elevation<-OurData$Elevation
# 
# if(exists("Latitude", OurData)) AnalysisDataFull$Latitude<-OurData$Latitude
# if(exists("Longitude", OurData)) AnalysisDataFull$Longitude<-OurData$Longitude
# if(exists("Easting", OurData)) AnalysisDataFull$Easting<-OurData$Easting
# if(exists("Northing", OurData)) AnalysisDataFull$Northing<-OurData$Northing
# 
# if(exists("Recording_Schedule", OurData)) AnalysisDataFull$Recording_Schedule<-OurData$Recording_Schedule
# if(exists("Mic_Type", OurData)) AnalysisDataFull$Mic_Type<-OurData$Mic_Type
# if(exists("file_duration_sec", OurData)) AnalysisDataFull$file_duration_sec<-OurData$file_duration_sec
# if(exists("file", OurData)) AnalysisDataFull$file<-OurData$file
# if(exists("audit_file", OurData)) AnalysisDataFull$audit_file<-OurData$audit_file
# if(exists("probeId", OurData)) AnalysisDataFull$probeId<-OurData$probeId
# if(exists("key", OurData)) AnalysisDataFull$key<-OurData$key
# 
# if(exists("Period", OurData)) AnalysisDataFull$Period<-OurData$Period
# if(exists("Round", OurData)) AnalysisDataFull$Round<-OurData$Round
# 
# if(exists("DateTime", OurData)) AnalysisDataFull$DateTime<-OurData$DateTime
# if(exists("Date", OurData)) AnalysisDataFull$Date<-OurData$Date
# if(exists("DayOfYear", OurData)) AnalysisDataFull$DayOfYear<-OurData$DayOfYear
# if(exists("NightOfYear", OurData)) AnalysisDataFull$NightOfYear<-OurData$NightOfYear
# if(exists("year", OurData)) AnalysisDataFull$year<-OurData$year
# if(exists("month", OurData)) AnalysisDataFull$month<-OurData$month
# if(exists("day", OurData)) AnalysisDataFull$day<-OurData$day
# if(exists("hour", OurData)) AnalysisDataFull$hour<-OurData$hour
# if(exists("minute", OurData)) AnalysisDataFull$minute<-OurData$minute
# 
# if(exists("sunrise_time", OurData)) AnalysisDataFull$sunrise_time<-OurData$sunrise_time
# if(exists("sunset_time", OurData)) AnalysisDataFull$sunset_time<-OurData$sunset_time
# if(exists("MinFromSunset", OurData)) AnalysisDataFull$MinFromSunset<-OurData$MinFromSunset
# if(exists("MinFromSunrise", OurData)) AnalysisDataFull$MinFromSunrise<-OurData$MinFromSunrise
# 
# if(exists("moonrise_time", OurData)) AnalysisDataFull$moonrise_time<-OurData$moonrise_time
# if(exists("moonset_time", OurData)) AnalysisDataFull$moonset_time<-OurData$moonset_time
# if(exists("MinFromMoonset", OurData)) AnalysisDataFull$MinFromMoonset<-OurData$MinFromMoonset
# if(exists("MinFromMoonrise", OurData)) AnalysisDataFull$MinFromMoonrise<-OurData$MinFromMoonrise
# if(exists("Illu", OurData)) AnalysisDataFull$Illu<-OurData$Illu
# if(exists("Moonup", OurData)) AnalysisDataFull$Moonup<-OurData$Moonup  
# if(exists("MoonAltitude", OurData)) AnalysisDataFull$MoonAltitude<-OurData$MoonAltitude  
# 
# # if(exists("civil_dawn", OurData)) AnalysisDataFull$civil_dawn<-OurData$civil_dawn
# # if(exists("civil_dusk", OurData)) AnalysisDataFull$civil_dusk<-OurData$civil_dusk
# # if(exists("MinFromCivilDawn", OurData)) AnalysisDataFull$MinFromCivilDawn<-OurData$MinFromCivilDawn
# # if(exists("MinFromCivilDusk", OurData)) AnalysisDataFull$MinFromCivilDusk<-OurData$MinFromCivilDusk
# # 
# # if(exists("naut_dawn", OurData)) AnalysisDataFull$naut_dawn<-OurData$naut_dawn
# # if(exists("naut_dusk", OurData)) AnalysisDataFull$naut_dusk<-OurData$naut_dusk
# # if(exists("MinFromNautDawn", OurData)) AnalysisDataFull$MinFromNautDawn<-OurData$MinFromNautDawn
# # if(exists("MinFromNautDusk", OurData)) AnalysisDataFull$MinFromNautDusk<-OurData$MinFromNautDusk
# # 
# # if(exists("astronomical_dawn", OurData)) AnalysisDataFull$astronomical_dawn<-OurData$astronomical_dawn
# # if(exists("astronomical_dusk", OurData)) AnalysisDataFull$astronomical_dusk<-OurData$astronomical_dusk
# # if(exists("MinFromAstroDawn", OurData)) AnalysisDataFull$MinFromAstroDawn<-OurData$MinFromAstroDawn
# # if(exists("MinFromAstroDusk", OurData)) AnalysisDataFull$MinFromAstroDusk<-OurData$MinFromAstroDusk
# 
# if(exists("flux", OurData)) AnalysisDataFull$flux<-OurData$flux
# if(exists("flux_sensitive", OurData)) AnalysisDataFull$flux_sensitive<-OurData$flux_sensitive
# if(exists("burst", OurData)) AnalysisDataFull$burst<-OurData$burst
# if(exists("level", OurData)) AnalysisDataFull$level<-OurData$level
# if(exists("level_absolute", OurData)) AnalysisDataFull$level_absolute<-OurData$level_absolute
# if(exists("click", OurData)) AnalysisDataFull$click<-OurData$click
# 
# if(exists("Total_alarms", OurData)) AnalysisDataFull$Total_alarms<-OurData$Total_alarms
# if(exists("Total_hits", OurData)) AnalysisDataFull$Total_hits<-OurData$Total_hits

# Data Check 4 (Final Check) --------------------------------------------------------------------------------

# DataCheck4<-AnalysisDataFull%>%
#     group_by(Sensor_Name)%>%
#     summarise(AnalysisDataFullTotals=sum(Total_hits))
# 
# DataCheck<-full_join(DataCheck,DataCheck4,by=c("location"="Sensor_Name"))
# DataCheck$diff3<-DataCheck$AnalysisDataFullTotals-DataCheck$OurDataTotals
# if(sum(DataCheck$diff3)!=0) stop("You are loosing or gaining detections")
# as.data.frame(DataCheck)

# Add 15 min time steps ---------------------------------------------------------------------------------

TimeStepOne = 17 # hour of first time step

getTimeStep<-function(baseTime, starthour, stepLength=15){
  require(lubridate)
  require(stringr)
  time_in_mins<-as.numeric(hour(baseTime))*60 + as.numeric(minute(baseTime))
  TimeStep<-time_in_mins%/%stepLength
  TimeStep<-TimeStep-(starthour*60)%/%stepLength
  TimeStep<-ifelse(TimeStep<0, TimeStep+((24*60)%/%stepLength), TimeStep)
  return(TimeStep)
}

getTimeStepStart<-function(baseTime, starthour, stepLength=15) {
  require(lubridate)
  require(stringr)
  time_in_mins<-as.numeric(hour(baseTime))*60 + as.numeric(minute(baseTime))
  TimeStep<-time_in_mins%/%stepLength
  TimeStepStart<-str_replace(paste(TimeStep%/%(60/stepLength), (TimeStep%%(60/stepLength))*stepLength, '00', sep=':'), ':0:', ':00:')
  return(TimeStepStart)
}

# Make time step
OurData$TimeStep<-getTimeStep(OurData$DateTime, TimeStepOne)
OurData$TimeStepStart<-getTimeStepStart(OurData$DateTime, TimeStepOne)
OurData$TimeStepN<-OurData$TimeStep-min(OurData$TimeStep)+1

OurData$SunsetTimeStep<-getTimeStep(OurData$sunset_time, TimeStepOne)
OurData$SunsetTimeStepN<-OurData$SunsetTimeStep-min(OurData$TimeStep)+1

OurData$SunriseTimeStep<-getTimeStep(OurData$sunrise_time, TimeStepOne)
OurData$SunriseTimeStepN<-OurData$SunriseTimeStep-min(OurData$TimeStep)+1

# Save AnalysisDataFull  ----------------------------------------------------------------------

sum(is.na(OurData$rating))

OurData$BinID<-paste(OurData$location,OurData$Date,OurData$TimeStepN)

OurDataPath<-'D:/CM,Inc/Dropbox (CMI)/CMI_Team/Analysis/Audit_Subset_Workflow/OurData_Subset_Audited.csv'
write.csv(OurData,file=OurDataPath,row.names=F)

# AnalysisDataPathRda<-file.path(ProjectDir,"R_output",paste("AnalysisData_",ProjectName,"_",Species,Today,".Rda",sep=""))
# saveRDS(AnalysisDataFull,file = AnalysisDataPathRda)
# 
# # AnalysisDataPathCSV<-file.path(ProjectDir,"R_output",paste("AnalysisData_",ProjectName,"_",Species,Today,".csv",sep=""))
# # write.csv(AnalysisDataFull,AnalysisDataPathCSV,row.names=F)
# 
# Totalend.time <- Sys.time()
# Totaltime.taken <- start.time - Totalend.time
# Totaltime.taken

######## JOHNSTONE STRAIT PURSE SEINE DATA CLEAN
###### Clean and configure raw data JS purse seine cruises for summer 2015
###### MAY 2016


rm(list = ls(all=TRUE)); #Remove all the objects in the memory


##########################################################
########## LIBRARIES ##########
library("plyr")
library("fossil")

################################

setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")


FishJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/IndividualFishDataJS2015.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
BridgeJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/BridgeLog_JS2015.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
OtolithJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/Otoliths_JS2015.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
# OtolithJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/TempOtoliths_JS2015.csv", 
#           stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))


##########################################################


########## INITIAL CLEAN ##########

#### Add year data
OtolithJS$Year <- 2015

#### Remove samples where otoliths not appropriate
OtolithJS <- OtolithJS[which(OtolithJS$Use == "Y"),]


## Edit universal number in otolith file to match fish file
temp <- strsplit(OtolithJS$UniversalNumber, "E")
OtolithJS$UniFishNumber <- unlist(temp)[2*(1:nrow(OtolithJS))]

FishJS <- plyr::rename(FishJS, c("universal.fish.number"="UniFishNumber"))


#### Trim each dataset before merging
head(FishJS)
FishJS_Trim <- FishJS[,c(3,4,8,11,24,26)]
head(FishJS_Trim)

head(BridgeJS)
BridgeJS_Trim <- BridgeJS[,c(1,2,10:13,19)]
head(BridgeJS_Trim)

head(OtolithJS)
OtolithJS_Trim <- OtolithJS[,c(1,5,9:23,27)]
# OtolithJS_Trim <- OtolithJS[,c(1,6,12:24,28)]
head(OtolithJS_Trim)


#### Merge datasets
temp <- merge(FishJS_Trim, BridgeJS_Trim, by="TOW_NUMBER")
DataJS2015 <- merge(temp, OtolithJS_Trim, by="UniFishNumber")


#### Rename variables
DataJS2015 <- plyr::rename(DataJS2015, c("TOW_NUMBER"="SetNumber",
	"Date.fished"="Date", "LENGTH"="ShipFL", "Stock.1"="Stock",
	"sub.area"="SubArea"))


#### Consolidate lat/long into single columns
DataJS2015$Lat <- (DataJS2015$START_LATITUDE_MINUTE/60) + DataJS2015$START_LATITUDE_DEGREE
DataJS2015$Long <- -1*((DataJS2015$START_LONGITUDE_MINUTE/60) + DataJS2015$START_LONGITUDE_DEGREE)


#### Add Julian date
DataJS2015$JulianDate <- strptime(DataJS2015$Date, "%d/%b/%Y")$yday + 1


#### Remove samples with probability <0.5
DataJS2015 <- DataJS2015[which(DataJS2015$Prob.1 > 0.49),]


table(DataJS2015$Stock)


#### Remove unnecessary columns and reorder
DataJS2015 <- DataJS2015[,c(1,2,3,32,13,12,30,31,5,4,18,20:29)]
# TempDataJS2015 <- DataJS2015[,c(1,2,30,4,14,15,18:27)]


##### Create and export file for use in BackCalc.R script which
#### supplies parameters for entryFL estimates below
write.csv(DataJS2015, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/JS2015_CF.csv")

# write.csv(TempDataJS2015, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/TempJS2015_CF.csv")


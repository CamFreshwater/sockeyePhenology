######## JOHNSTONE STRAIT PURSE SEINE DATA CLEAN
###### Clean and configure raw data JS purse seine cruises for summer 2016
###### MAY 2017


rm(list = ls(all=TRUE)); #Remove all the objects in the memory


##########################################################
########## LIBRARIES ##########
library("plyr")
library("fossil")

################################

setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")


FishJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2016/DIPS2016_CamSelect.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
BridgeJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2016/BridgeLog_JS2016.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
OtolithJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2016/Otoliths_JS2016.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))



##########################################################


########## INITIAL CLEAN ##########

#### Add year data
OtolithJS$Year <- 2016

#### Remove samples where otoliths not appropriate
OtolithJS <- OtolithJS[which(OtolithJS$Use == "Y"),]


OtolithJS <- plyr::rename(OtolithJS, c("UniversalNumber"="UniFishNumber"))
FishJS <- plyr::rename(FishJS, c("universal.fish.number"="UniFishNumber"))



#### Trim each dataset before merging
str(FishJS)
FishJS_Trim <- FishJS[,c(2,3,5,7,16,18)]
head(FishJS_Trim)

str(BridgeJS)
BridgeJS_Trim <- BridgeJS[,c(1,2,19,20,22,23)]
head(BridgeJS_Trim)

str(OtolithJS)
OtolithJS_Trim <- OtolithJS[,c(1,6,12:24)]
head(OtolithJS_Trim)


#### Merge datasets
temp <- merge(FishJS_Trim, BridgeJS_Trim, by="TOW_NUMBER")
DataJS2016 <- merge(temp, OtolithJS_Trim, by="UniFishNumber")


#### Rename variables
DataJS2016 <- plyr::rename(DataJS2016, c("TOW_NUMBER"="SetNumber",
	"LENGTH"="ShipFL", "Stock.1"="Stock"))
head(DataJS2016)


#### Consolidate lat/long into single columns
DataJS2016$Lat <- (DataJS2016$START_LATITUDE_MINUTE/60) + DataJS2016$START_LATITUDE_DEGREE
DataJS2016$Long <- -1*((DataJS2016$START_LONGITUDE_MINUTE/60) + DataJS2016$START_LONGITUDE_DEGREE)


#### Add Julian date
DataJS2016$JulianDate <- strptime(DataJS2016$Date, "%d-%b-%Y")$yday + 1


#### Remove samples with probability <0.5
DataJS2016 <- DataJS2016[which(DataJS2016$Prob.1 > 0.49),]


table(DataJS2016$Stock)


#### Remove unnecessary columns and reorder
DataJS2016 <- DataJS2016[,c(1,2,28,26,27,4,5,12,16:24)]


#### TEMPORARY - remove sets 31/71 until lat/long can be estimated
DataJS2016 <- DataJS2016[complete.cases(DataJS2016$Lat),]


##### Create and export file for use in BackCalc.R script which
#### supplies parameters for entryFL estimates below
write.csv(DataJS2016, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2016/JS2016_CF.csv")








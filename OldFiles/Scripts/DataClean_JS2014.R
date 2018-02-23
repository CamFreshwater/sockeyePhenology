######## JOHNSTONE STRAIT PURSE SEINE DATA CLEAN
###### Clean and configure raw data JS purse seine cruises starting summer of 2014
###### MAY 2015


rm(list = ls(all=TRUE)); #Remove all the objects in the memory



setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")


GeneticsJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2014/IndividualStockIDs_JS2014.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
FishJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2014/IndividualFishData_JS2014.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
BridgeJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2014/BridgeLog_JS2014.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
OtolithJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2014/Otoliths_JS2014.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))


##########################################################
########## LIBRARIES ##########
library("plyr")
library("fossil")
##########################################################


########## INITIAL CLEAN ##########


#### Trim each dataset before merging
head(GeneticsJS)
GeneticsJS_Trim <- GeneticsJS[,c(1,2,11:13)] 
GeneticsJS_Trim <- rename(GeneticsJS_Trim, c("VIAL.." = "DNA"))

head(FishJS)
FishJS_Trim <- FishJS[,c(3,12,28,33,37)]

head(BridgeJS)
BridgeJS_Trim <- BridgeJS[,c(1,2,22:25,39)]

head(OtolithJS)
OtolithJS_Trim <- OtolithJS[,c(1,5,9:18)]


#### Merge datasets and rename variables
temp <- merge(FishJS_Trim, BridgeJS_Trim, by="TOW_NUMBER")
ShipDataJS2014 <- merge(temp, GeneticsJS_Trim, by="DNA")

#### Rename variables
ShipDataJS2014 <- rename(ShipDataJS2014, c("OTHER_WEIGHT" = "LabWeight"))
ShipDataJS2014 <- rename(ShipDataJS2014, c("LENGTH" = "ShipFL"))
ShipDataJS2014 <- rename(ShipDataJS2014, c("EVENT_DATE" = "Date"))
ShipDataJS2014 <- rename(ShipDataJS2014, c("AREA_CODE" = "AreaCode"))
ShipDataJS2014 <- rename(ShipDataJS2014, c("TOW_NUMBER" = "Tow"))
ShipDataJS2014 <- rename(ShipDataJS2014, c("FISH_ID" = "FishNumber"))
ShipDataJS2014 <- rename(ShipDataJS2014, c("Stock.1" = "Stock"))

DataJS2014 <- merge(ShipDataJS2014, OtolithJS_Trim, by="FishNumber")


#### Consolidate lat/long into single columns
DataJS2014$Lat <- (DataJS2014$START_LATITUDE_MINUTE/60) + DataJS2014$START_LATITUDE_DEGREE
DataJS2014$Long <- -1*((DataJS2014$START_LONGITUDE_MINUTE/60) + DataJS2014$START_LONGITUDE_DEGREE)

#### Add julian dates
DataJS2014$JulianDate <- strptime(DataJS2014$Date, "%m/%d/%Y")$yday+1


#### Remove extra columns and reorder
DataJS2014 <- DataJS2014[,c(1,6,29,27,28,11,3,4,5,2,13:15,17:26)]
#### Remove adult fish (i.e. shipFL >250)
DataJS2014 <- DataJS2014[which(DataJS2014$ShipFL < 250),]
#### Remove samples where otoliths not appropriate
DataJS2014 <- DataJS2014[which(DataJS2014$Use == "Y"),]
#### Remove samples where stock ID probability <0.50
DataJS2014 <- DataJS2014[which(DataJS2014$Prob.1 > 0.49),]


#### Convert variables and exploratory analyses
DataJS2014$Stock <- as.factor(DataJS2014$Stock)
table(DataJS2014$Stock) 
# Retain samples from populations with at least 8 individuals
DataJS2014 <- DataJS2014[which(DataJS2014$Stock == "Birkenhead" | DataJS2014$Stock == "Chilko" | 
	DataJS2014$Stock == "DollyVarden_Cr" | DataJS2014$Stock == "Nadina" | DataJS2014$Stock == "Pitt"),]
DataJS2014$AreaCode <- as.factor(DataJS2014$AreaCode)
table(DataJS2014$AreaCode) 
DataJS2014$Region.1 <- as.factor(DataJS2014$Region.1)
table(DataJS2014$Region.1) 
#############


##### Create and export file for use in BackCalc.R script which
#### supplies parameters for entryFL estimates below
write.csv(DataJS2014, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2014/JS2014_CF.csv")

################################################################
################################################################
################################################################


##### FOLLOWING DATA NOW CALCULATED IN DataCombineJS.R script



#### Back calculate estimates of fork length at ocean entry average of scale-proportional and body-proportional
### hypotheses from Francis 1990
### Use regressions from BackCalc.R that have all sockeye in them


##### Convert oto metrics to mm so comparable with FL

# DataJS2014$OtoRadius <- (DataJS2014$Total_Radius)/1000
# DataJS2014$MECWidth <- (DataJS2014$MEC_Radius)/1000


# # ##### Scale Proportional Estimate 

# ### 1) Calculate otolith width based on fork length
# OW.lm <- lm(OtoRadius ~ ShipFL, data=DataJS2014) 
# summary(OW.lm)
# ### 2) Calculate mean OW for each fish based on regression
# # Regression: OW = 0.002693FL + 0.3699
# OW.PopCapture <- (0.0026974*DataJS2014$ShipFL) + 0.3672
# ### 3) Calculate ratio of observed to expected
# OW.Ratio <- DataJS2014$OtoRadius/OW.PopCapture
# ### 4) Adjust observed MEC width by ratio to calculate expected
# OW.ExpMEC <- DataJS2014$MECWidth/OW.Ratio
# ### 5) Estimate entry length based on proportion
# DataJS2014$OW.EntryFL <- (OW.ExpMEC - 0.3672)/0.0026974


# # ##### Length Proportional Estimate 

# ### 1) Calculate FL based on OW
# FL.lm <- lm(ShipFL ~ OtoRadius, data=DataJS2014) 
# summary(FL.lm)
# ### 2) Calculate mean FL for each fish based on regression
# # Regression: FL = 290.159W - 79.315
# FL.PopCapture <- (293.076*DataJS2014$OtoRadius) - 80.799
# ### 3) Calculate ratio of observed to expected
# FL.Ratio <- DataJS2014$ShipFL/FL.PopCapture
# ### 4) Calculate mean width for individual with observed MEC width
# FL.MECwidth <- (293.076*DataJS2014$MECWidth) - 80.799
# ### 5) Adjust by ratio
# DataJS2014$FL.EntryFL <- FL.Ratio * FL.MECwidth


# # #### Take mean of both methods to reduce error
# DataJS2014$EntryFL <- (DataJS2014$FL.EntryFL + DataJS2014$OW.EntryFL)/2
# #######################################




# #### Add migration distance by converting difference in lat-long between river mouth 
# ### and capture location; since all stocks with sufficient sample size are FR, 
# ### only one start lat-long  (49.083 N, -123.2 W) 
# DataJS2014$MigDist <- deg.dist(-123.2,49.083,DataJS2014$Long,DataJS2014$Lat)
     
# ### Estimate migration rate in km/day
# DataJS2014$MigRate <- DataJS2014$MigDist/DataJS2014$Total_Count

# ### Estimate migration rate in body lengths (at entry) per second
# # Convert distance from km to body lengths
# MigDistBL <- (DataJS2014$MigDist*1000000)/DataJS2014$EntryFL
# # Convert days to seconds and calculate rate
# DataJS2014$MigRateBLS <- MigDistBL/(DataJS2014$Total_Count*86400)


# #### Add mean growth rate
# DataJS2014$BodyGrowth <- (DataJS2014$ShipFL - DataJS2014$EntryFL)/DataJS2014$Total_Count
# DataJS2014$OtoGrowth <- DataJS2014$MarineGrowth/DataJS2014$Total_Count


# #### Add entry date 
# DataJS2014$EntryDate <- (DataJS2014$JulianDate - DataJS2014$Total_Count)


# write.csv(DataJS2014, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/2014/JS2014_CF.csv")

# #### Trim and clean .csv to send to Julia with growth data
# CropJS2014 <- DataJS2014[,c(1:11,14,15:22,28,34,32,33)]
# write.csv(CropJS2014, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/2014/JS2014_ForDFO.csv")




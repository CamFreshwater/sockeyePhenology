######## JOHNSTONE STRAIT PURSE SEINE COMBINE AND CLEAN
###### Combine 2014 and 2015 JS data
###### MAY 2016


rm(list = ls(all=TRUE)); #Remove all the objects in the memory



setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")


JS2014 <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2014/JS2014_CF.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
JS2015 <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2015/JS2015_CF.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
JS2016 <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Data2016/JS2016_CF.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))


require(plyr)
require(fossil)



################# DATA CLEAN #################


#### Add year variables
JS2014$Year <- 2014


#### Trim to ship and otolith data only 
JS2014 <- JS2014[,c(2,4,25,5,6,12,8,9,17,15,16,18:23)]
JS2015 <- JS2015[,c(2,5,6,8:10,3,11,14:22)]
JS2016 <- JS2016[,c(2,4,9,5,6,8,3,7,10:18)]


#### Rename variables
JS2014 <- plyr::rename(JS2014, c("MEC_Radius"="MECRadius",
	"Total_Radius"="TotalRadius", "Total_Count"="TotalCount", "MarineGrowth"="MarGrowth"))
JS2015 <- plyr::rename(JS2015, c("UniFishNumber"="FishNumber",
	"SetNumber"="Tow"))
JS2016 <- plyr::rename(JS2016, c("UniFishNumber"="FishNumber",
	"SetNumber"="Tow"))

#### Combine datasheets
DataJS <- rbind(JS2014, JS2015, JS2016)


#### Add CU level based on spawning pop
DataJS$CU <- "Quesnel-S"
DataJS$CU[DataJS$Stock == "Birkenhead"] <- "Lillooet-Harrison-L"
DataJS$CU[DataJS$Stock == "Chilko" | DataJS$Stock == "Chilko_south"] <- "Chilko-S/Chilko-ES"
DataJS$CU[DataJS$Stock == "DollyVarden_Cr"] <- "Chilliwack-ES"
DataJS$CU[DataJS$Stock == "L_Adams" | DataJS$Stock == "L_Shuswap" | 
	DataJS$Stock == "MiddleShuswap"] <- "Shuswap-L"
DataJS$CU[DataJS$Stock == "Seymour"] <- "Shuswap-ES"
DataJS$CU[DataJS$Stock == "Nadina"] <- "Nadina-Francois-ES"
DataJS$CU[DataJS$Stock == "Pitt"] <- "Pitt-ES"
DataJS$CU[DataJS$Stock == "Stellako"] <- "Francois-Fraser-S"


#### Add migration distance covariate based on nursery lake
DataJS$RiverDist <- 736
DataJS$RiverDist[DataJS$Stock == "Birkenhead"] <- 185.5
DataJS$RiverDist[DataJS$Stock == "Chilko" | DataJS$Stock == "Chilko_south"] <- 661
DataJS$RiverDist[DataJS$Stock == "DollyVarden_Cr"] <- 152.07
DataJS$RiverDist[DataJS$Stock == "L_Adams" | DataJS$Stock == "L_Shuswap" | 
	DataJS$Stock == "MiddleShuswap"] <- 524
DataJS$RiverDist[DataJS$Stock == "Seymour"] <- 517
DataJS$RiverDist[DataJS$Stock == "Nadina"] <- 1076
DataJS$RiverDist[DataJS$Stock == "Pitt"] <- 55
DataJS$RiverDist[DataJS$Stock == "Stellako"] <- 986


#### Add latitude distance covariate based on nursery lake
DataJS$LakeLat <- 52.51667
DataJS$LakeLat[DataJS$Stock == "Birkenhead"] <- 49.9
DataJS$LakeLat[DataJS$Stock == "Chilko" | DataJS$Stock == "Chilko_south"] <- 51.26667
DataJS$LakeLat[DataJS$Stock == "DollyVarden_Cr"] <- 49.05
DataJS$LakeLat[DataJS$Stock == "L_Adams" | DataJS$Stock == "L_Shuswap" | 
	DataJS$Stock == "MiddleShuswap"] <- 50.89
DataJS$LakeLat[DataJS$Stock == "Seymour"] <- 50.98
DataJS$LakeLat[DataJS$Stock == "Nadina"] <- 53.95833
DataJS$LakeLat[DataJS$Stock == "Pitt"] <- 49.43333
DataJS$LakeLat[DataJS$Stock == "Stellako"] <- 54.05833



# DataJS <- DataJS[,c(1:6,18,19,7:17)]


#### Export file for BackCalculation
# write.csv(DataJS, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/FullJS.csv")


DataJS$Stock <- as.factor(DataJS$Stock)
DataJS$CU <- as.factor(DataJS$CU)
DataJS$Year <- as.factor(DataJS$Year)
DataJS$Age <- as.factor(DataJS$Age)


#### Add variables
DataJS$EntryDate <- round(DataJS$JulianDate - DataJS$TotalCount)
DataJS$OtoGrowth <- DataJS$MarGrowth / DataJS$TotalCount


#######################################
######### Calculate entry FL

##### Convert oto metrics to mm so comparable with FL
DataJS$OtoRadius <- (DataJS$TotalRadius)/1000
DataJS$MECWidth <- (DataJS$MECRadius)/1000



########### USES GLOBAL VALUES ###########

#### CURRENT REGRESSION PARAMETERS FOR BACK CALC IN TOTAL (don't work for age-0)
# Regression (UPDATED JAN 24 2017): OW = 0.002696FL + 0.3682
# Regression (UPDATED JAN 24 2017): FL = 285.103W - 76.195

## 1) Calculate otolith width based on fork length
# OW.lm <- lm(ShipFL ~ OtoRadius, data=DataJS) 
# summary(OW.lm)

# ### 2) Calculate mean OW for each fish based on regression
# OW.PopCapture <- (0.002696*DataJS$ShipFL) + 0.3682
# ### 3) Calculate ratio of observed to expected
# OW.Ratio <- DataJS$OtoRadius/OW.PopCapture
# ### 4) Adjust observed MEC width by ratio to calculate expected
# OW.ExpMEC <- DataJS$MECWidth/OW.Ratio
# ### 5) Estimate entry length based on proportion
# DataJS$OW.EntryFL <- (OW.ExpMEC - 0.3682)/0.002696

# # # ##### Length Proportional Estimate 

# ### 1) Calculate FL based on OW
# # FL.lm <- lm(ShipFL ~ OtoRadius, data=DataJS) 
# # summary(FL.lm)
# ### 2) Calculate mean FL for each fish based on regression
# FL.PopCapture <- (285.103*DataJS$OtoRadius) - 76.195
# ### 3) Calculate ratio of observed to expected
# FL.Ratio <- DataJS$ShipFL/FL.PopCapture
# ### 4) Calculate mean width for individual with observed MEC width
# FL.MECwidth <- (285.103*DataJS$MECWidth) - 76.195
# ### 5) Adjust by ratio
# DataJS$FL.EntryFL <- FL.Ratio * FL.MECwidth

# # # #### Take mean of both methods to reduce error
# DataJS$EntryFL2 <- (DataJS$FL.EntryFL + DataJS$OW.EntryFL)/2
###############################################


########### ESTIMATE ENTRY SIZE WITH BIOLOGICAL INTERCEPT ###########
L0 = 19.76 ## average alevin length from intermediate  treatments (Beacham and Murray 1988)
Lcpt = DataJS$ShipFL
Ri = DataJS$MECWidth
# R0 = mean(DataJS$InRadius)/1000 #0.226 (from temporary JS2015 DataJSaset)
R0 = 0.226
## average hatch check radius
Rcpt = DataJS$OtoRadius


### Plot of otolith size vs. ship size w/ proposed intercept 
plot(OtoRadius ~ ShipFL, data=DataJS, xlim=c(0,155), ylim=c(0,0.9))
points(19.76, 0.226, pch=16, col="red")
abline(lm(OtoRadius ~ ShipFL, data=DataJS))


### Calculate entry size
DataJS$EntryFL <- Lcpt + ((Ri - Rcpt) * (Lcpt - L0) * (Rcpt - R0)^-1)


hist(DataJS$EntryFL)
plot(EntryFL ~ ShipFL, data=DataJS)


#### What do predicted SEs look like?
reg <- lm(Lcpt ~ Rcpt + DataJS$Age)
summary(reg)

z <- predict(reg, newdata=data.frame(Age=DataJS$Age, Rcpt=Rcpt), se.fit=TRUE)
z$se.fit

#### Use NLS to fit above model and get SEs of predictions
fit <- nls(Lcpt + ((Ri - Rcpt) * (Lcpt - L0) * (Rcpt - R0)^-1))


#######################################
#######################################

### Add growth data
DataJS$BodyGrowth <- (DataJS$ShipFL - DataJS$EntryFL) / DataJS$TotalCount
hist(DataJS$BodyGrowth)
### Scaled body growth by initial size 
### %growth per day
DataJS$ScaledGrowth <- (DataJS$BodyGrowth/DataJS$EntryFL)*100


## Add migration distance by converting difference in lat-long between river mouth 
# and capture location; since all stocks with sufficient sample size are FR, 
# only one start lat-long  (49.083 N, -123.2 W) 
DataJS$MigDist <- deg.dist(-123.2,49.083,DataJS$Long,DataJS$Lat)
     
### Estimate migration rate in km/day
DataJS$MigRate <- DataJS$MigDist/DataJS$TotalCount

### Estimate migration rate in body lengths (at entry) per second
# Convert distance from km to body lengths
MigDistBL <- (DataJS$MigDist*1000000)/DataJS$EntryFL
# Convert days to seconds and calculate rate
DataJS$MigRateBLS <- MigDistBL/(DataJS$TotalCount*86400)


### Add midpoint day
DataJS$MidDay <- (DataJS$TotalCount/2) + DataJS$EntryDate


#######################################################
#######################################################
write.csv(DataJS, file="/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv")


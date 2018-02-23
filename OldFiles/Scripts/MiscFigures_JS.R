##################################################################
##################################################################
######## Johnstone Strait - Figures 
######## January 2017
##################################################################
##################################################################


rm(list = ls(all=TRUE)); #Remove all the objects in the memory


setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")


JSdat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
### Import lake dataset file from DFA analysis for labeling nursery lakes
LakeDat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/LakeData.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))



JSdat <- JSdat[,-1]

####################################################################
############################ LIBRARIES ###############################
####################################################################

library(mapdata)
library(maps)
library(maptools)
library(PBSmapping)
library(scales)
library(HDMD) #pairwise.mahalanobis function
library(cluster) #agnes function
library(vegan)
library(ecodist) #non pairwise mahalanobis function
library(plyr)
data(nepacLL)


######################################################################

####################################################################
############################# CLEAN ################################
####################################################################

JSdat$Year <- as.factor(JSdat$Year)
JSdat$Age <- as.factor(JSdat$Age)


#### Lake data

### Trim to relevant CUs, keeping only primary lakes
temp <- LakeDat[LakeDat$CU == "Chilko" | LakeDat$CU == "Francois-Fraser" | LakeDat$CU == "Chilliwack" |
	 LakeDat$CU == "Lillooet-Harrison" | LakeDat$CU == "Nadina-Francois" | LakeDat$CU == "Pitt" | 
	 LakeDat$CU == "Quesnel" | LakeDat$CU == "Shuswap-L",]
LakeDat <- temp[temp$Imp=="1",]


LakeDat$CU <- as.factor(LakeDat$CU)
LakeDat$RearingLake <- as.factor(LakeDat$RearingLake)

#### Merge lat and long measurements
LakeDat$Lat <- LakeDat$Lat1 + (LakeDat$Lat2/60)
LakeDat$Long <- LakeDat$Long1 + (LakeDat$Long2/60)


################# FIGURE 1: Watershed map #################

### Read in shape file
# lakes=readShapePoly("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/Environmental/Shapefiles/lakes_bc_dd/lakes_bc_dd.shp")
coast=readShapeSpatial("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/Environmental/Shapefiles/coastline_bc_dd/coastline_bc_dd.shp")
rivers=readShapePoly("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/Environmental/Shapefiles/BCRivers/FWRVRSPL_polygon.shp")
lakes=readShapePoly("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/Environmental/Shapefiles/FIS_LK_SVW/FIS_LK_SVW_polygon.shp")



###### Clean shape files for plotting by removing extra lakes and rivers

### Remove lake data without a gazetted name
temp <- lakes[complete.cases(lakes@data[,6]),]

### Exclude all lakes not included in analysis
temp.lak <- temp[c(temp@data$GZTTDNM=="SHUSWAP LAKE" | 
	# temp@data$GZTTDNM=="BOWRON LAKE" | temp@data$GZTTDNM=="STUART LAKE" |
	temp@data$GZTTDNM=="CHILKO LAKE" | temp@data$GZTTDNM=="FRASER LAKE" |
	temp@data$GZTTDNM=="FRANCOIS LAKE" | temp@data$GZTTDNM=="HARRISON LAKE" |
	temp@data$GZTTDNM=="KAMLOOPS LAKE" | temp@data$GZTTDNM=="LILLOOET LAKE" |
	# temp@data$GZTTDNM=="NORTH BARRIERE LAKE" | temp@data$GZTTDNM=="ANDERSON LAKE" | temp@data$GZTTDNM=="SETON LAKE" |
	temp@data$GZTTDNM=="PITT LAKE" |
	temp@data$GZTTDNM=="QUESNEL LAKE" | temp@data$GZTTDNM=="ADAMS LAKE" |
	temp@data$GZTTDNM=="LITTLE SHUSWAP LAKE" | temp@data$GZTTDNM=="MARA LAKE" | 
	temp@data$GZTTDNM=="MABEL LAKE" | temp@data$GZTTDNM=="NADINA LAKE" |
	# temp@data$GZTTDNM=="TAKLA LAKE" | temp@data$GZTTDNM=="TREMBLEUR LAKE" | temp@data$GZTTDNM=="STUART LAKE"|
	temp@data$GZTTDNM=="CHILLIWACK LAKE" ),]
################


### Remove river data without a gazetted name
temp <- rivers[complete.cases(rivers@data[,13]),]

### Exclude all rivers not included in analysis
temp.riv <- temp[c(temp@data$GNSNM1=="Pitt River" | temp@data$GNSNM1=="Harrison River" |
	temp@data$GNSNM1=="Lillooet River" | temp@data$GNSNM1=="Thompson River" |
	# temp@data$GNSNM1=="Gates Creek" | temp@data$GNSNM1=="Seton River" | 
	temp@data$GNSNM1=="North Thompson River" | temp@data$GNSNM1=="South Thompson River" |
	temp@data$GNSNM1=="Shuswap River" | temp@data$GNSNM1=="Adams River" |
	temp@data$GNSNM1=="Quesnel River" | 
	# temp@data$GNSNM1=="Bowron River" |
	temp@data$GNSNM1=="Chilko River" | temp@data$GNSNM1=="Chilcotin River" |
	temp@data$GNSNM1=="Stellako River" | temp@data$GNSNM1=="Nechako River" |
	# temp@data$GNSNM1=="Stuart River" | temp@data$GNSNM1=="Tachie River" |temp@data$GNSNM1=="Middle River" | 
	temp@data$GNSNM1=="Fraser River" |
	temp@data$GNSNM1=="Chilliwack River" | temp@data$GNSNM1=="Chilliwack Creek" |
	temp@data$GNSNM1=="Nadina River"
	),]




setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Figures")



## set up plot space for main map [1] + inset [2]
# mat<-matrix(1, 50, 50)
# mat[24:49, 1:25]<-2
# mat[1:50, 26:50]<-3

mat<-matrix(1, 100, 100)
mat[65:100, 1:100]<-2

range(JSdat$Lat)
range(JSdat$Long)

png("StudyMap.png", width=4, height=6.25, units="in", res=700)

layout(mat)


#### Upper map showing watershed with nursery lakes
par(mar=c(1,0.2,0.2,0)+0.1, oma=c(0,0,0,0)+0.1,  cex.axis=1.1)
plotMap(nepacLL, xlim=c(-128,-118), ylim=c(48,55), 
	plt=NULL, col="lightgrey", bg="white", 
	cex=1, axes=FALSE, ylab="", xlab="")
plot(temp.lak, xlim=c(-126,-118.25), ylim=c(48.25,56), col="black", add=TRUE)
plot(temp.riv, xlim=c(-125.5,-118.25), ylim=c(48.25,55.25), lwd=0.75, col="white", 
	add=TRUE)
axis(1, tick=T, at=c(-130,-127,-123,-119,-116), lwd=2, lwd.ticks=2, cex.axis=1.25,
  labels=c("",expression(paste("127 ",degree,"W")),
    expression(paste("123 ",degree,"W")),expression(paste("119 ",degree,"W")),
    ""), mgp=c(2.4, 0.3, 0), tck=-0.02) 
axis(2, tick=T, at=c(47,50,54,57), lwd=2, lwd.ticks=2, cex.axis=1.25,
  labels=c("",expression(paste("50 ",degree,"N")),expression(paste("54 ",degree,"N")),""),
  mgp=c(2.4, 0.3, 0), tck=-0.02) 
axis(3, tick=T, at=c(-130,-117), lwd=2, labels=FALSE)
axis(4, tick=T, at=c(47,57), lwd=2, labels=FALSE)
polygon(x=c(-126.22,-125,-125,-126.22), y=c(50,50,50.5,50.5), lwd=1.45, border="red")


### Extract lat/long from lake dataset to use for plotting ID numbers
long <- -1*LakeDat$Long
lat <- LakeDat$Lat

### Add numbers for different lakes
text("1", x=long[4]-0.4, y=lat[4]+0.25, cex=1.4) #Francois
text("2", x=long[2]-0.2, y=lat[2]+0.25, cex=1.4) #Fraser
text("3", x=long[6]-0.25, y=lat[6]+0.2, cex=1.4) #Quesnel
text("4", x=long[7]-0.3, y=lat[7]+0.3, cex=1.4) #early Shuswap
text("5", x=long[7]+0.4, y=lat[7]+0.1, cex=1.4) #main Shuswap
text("6", x=long[1]+0.2, y=lat[1]+0.15, cex=1.4) #Chilko
text("7", x=long[3]+0.35, y=lat[3]-0.1, cex=1.4) #Lillooet
text("8", x=long[5]-0.3, y=lat[5]+0.15, cex=1.4) #Pitt
text("9", x=long[8]+0.2, y=lat[8]+0.15, cex=1.4) #Chilliwack


#### Lower map showing juvenile sampling locations
# par(mar=c(0,.2,0,.2)+0.1, oma=c(0,0,0,0)+0.1,  cex.axis=1.2)
plotMap(nepacLL, xlim=c(-126.33,-124.87), ylim=c(50,50.5), 
	plt=NULL, bg="white", col="lightgrey", 
	cex=1, axes=FALSE, ylab="", xlab="")
axis(1, tick=T, at=c(-127,-126,-125.2,-124.5), lwd=2, lwd.ticks=2, cex.axis=1.25,
  labels=c("",expression(paste("126 ",degree,"W")),expression(paste("125.2 ",degree,"W")),
    ""), mgp=c(2.4, 0.3, 0), tck=-0.04) 
axis(2, tick=T, at=c(49,50.1,50.4,51), lwd=2, lwd.ticks=2, cex.axis=1.25,
  labels=c("",expression(paste("50.1 ",degree,"N")), expression(paste("50.4 ",degree,"N")),""),
  mgp=c(2.4, 0.3, 0), tck=-0.04) 
axis(3, tick=T, at=c(-127,-124), lwd=2, labels=FALSE)
axis(4, tick=T, at=c(49,51), lwd=2, labels=FALSE)


### Add points from surveys
points(x=JSdat$Long, y=JSdat$Lat, bg=scales::alpha("red", 0.6), pch=21)


### Scale bar
map.scale(-126.15, 50.07, relwidth=0.14, metric=TRUE, ratio=FALSE, 
	cex=0.9)


dev.off()

######################################################################
######################################################################
######################################################################



####### FIGURE 2 - CLUSTER ANALYSIS ####### 

### Use pairwise Mahalanobis which allows groups to be easily recognized
## and accounts for original correlation structure of the variables
year = unique(JSdat$Year)


pdf("MigTraitCluster.pdf", width=5, height=2.5)

par(mfrow=c(1,3), mar=c(8,1.25,1.5,1)+0.1, oma=c(1,3,0,0))

for(i in 1:length(year)){
	temp <- subset(JSdat, JSdat$Year == year[i])
	
	cu = matrix(temp$CU)
	cu = t(cu[,1])

	mig_traits = c("TotalCount","BodyGrowth","EntryFL","EntryDate") 
	### NOTE: cannot include midday and entry date and calculate mahalanobis
	### distance; too strongly correlated; AC higher for entrydate so 
	### currently using those
	mig_traits = as.matrix(temp[,mig_traits])

	mahala_sq = pairwise.mahalanobis(x=mig_traits, grouping=cu)
	names = rownames(mahala_sq$means) #extract labels

	mahala = sqrt(mahala_sq$distance)
	rownames(mahala) = names
	colnames(mahala) = names

	cluster = agnes(mahala,diss=TRUE,keep.diss=FALSE,method="complete") #hierarchical clustering
	# plot(cluster, which.plots=2, main="", ylim=c(0,40)) #plot dendrogram
	plot(as.dendrogram(cluster), ylim=c(0,40))
	title(year[i])

	## store and calculate PERMANOVAs (note these use Euclidean
	## distances, not Mahalanobis) 
	# m1 <- adonis(mig_traits ~ as.factor(cu), permutations=999, method="euclidean")
	# assign(paste("mod", year[i], sep="."), m1)
	# print(m1) 
}

mtext(side=2, text=expression(paste("Mahalanobis Distance")), las=3, outer=TRUE,
 line=1)

dev.off()



#### Alternative cluster analysis that calculates group level means and variances,
### then groups CUs
### DOESN'T WORK BECAUSE EITHER UNSTABLE OR HCLUST CAN'T RECEIVE OUTPUT
year = unique(JSdat$Year)


for(i in 1:length(year)){
	temp <- subset(JSdat, JSdat$Year == year[i])
	
	cu = matrix(temp$CU)
	cu = t(cu[,1])

	traits.dat <- ddply(temp, .(as.factor(CU)), summarize, meanCount=mean(TotalCount), varCount=var(TotalCount), meanBodyGrowth=mean(BodyGrowth), 
		varBodyGrowth=var(BodyGrowth), meanEntryFL=mean(EntryFL), varEntryFL=var(EntryFL), meanEntryDate=mean(EntryDate), varEntryDate=var(EntryDate))
	traits.mat <- as.matrix(traits.dat[,-1])

	# # mig.dist <- distance(traits.mat, method="mahalanobis")
	# mig.dist <- mahalanobis(traits.mat, center=colMeans(traits.mat), cov(traits.mat), tol=9e-30)
	# mig.mahala <- hclust(mig.dist, method="average")
	# plot(mig.mahala, main="")
	
	# title(year[i])
	mahala_sq = pairwise.mahalanobis(x=traits.mat)

	mahala = sqrt(mahala_sq$distance)
	rownames(mahala) = cu
	colnames(mahala) = cu

	cluster = agnes(mahala,diss=TRUE,keep.diss=FALSE,method="complete") #hierarchical clustering
	# plot(cluster, which.plots=2, main="", ylim=c(0,40)) #plot dendrogram
	plot(as.dendrogram(cluster), ylim=c(0,40))
	title(year[i])
}

	temp <- subset(JSdat, JSdat$Year == year[3])
	
	cu = matrix(temp$CU)
	cu = t(cu[,1])

	traits.dat <- ddply(temp, .(as.factor(CU)), summarize, meanCount=mean(TotalCount), meanBodyGrowth=mean(BodyGrowth), 
		meanEntryFL=mean(EntryFL), meanEntryDate=mean(EntryDate))
	traits.mat <- as.matrix(traits.dat[,-1])

	mig.dist <- distance(traits.mat, method="mahalanobis")
	# mig.dist <- mahalanobis(traits.mat, center=colMeans(traits.mat), cov(traits.mat), tol=9e-30)
	mig.mahala <- hclust(mig.dist, method="average")
	plot(mig.mahala, main="")
	
	# # title(year[1])
	# mahala_sq = pairwise.mahalanobis(x=traits.mat)

	# mahala = sqrt(mahala_sq$distance)
	# rownames(mahala) = cu
	# colnames(mahala) = cu

	# cluster = agnes(mahala,diss=TRUE,keep.diss=FALSE,method="complete") #hierarchical clustering
	# # plot(cluster, which.plots=2, main="", ylim=c(0,40)) #plot dendrogram
	# plot(as.dendrogram(cluster), ylim=c(0,40))
	# title(year[1])


######################################################################
######################################################################
######################################################################



############# FIGURE 3/4 - KERNEL DENSITY PLOTS #############


## Make palette dataframe of colors for plotting based on region
## Gradient colors
# palette <- data.frame(CU=unique(d$CU), col=c("#d73027","#f46d43","#fdae61",
# 	"#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850"))
## Discrete colors
temp <- unique(JSdat[c("CU", "cu_id")])
temp <- temp[order(temp$cu_id),] # reorder based on CU id
palette <- data.frame(CU=temp$CU, cu_id=temp$cu_id, col=c("#a6cee3","#1f78b4","#b2df8a",
	"#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"))
palette$col<-as.character(palette$col)



####### WRITE FUNCTION #######

setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Figures")


### Loop to subset dataframe based on factor levels and 
## plot multiple density plots
pdf("CountKernel.pdf", width=6.65, height=6.65)

temp <- JSdat[,c("TotalCount","CU","Year")]

par(mfrow=c(3,1), mar=c(1.75,1.25,1.5,1)+0.1, oma=c(3,3,0,0), mgp=c(3,0.5,0), las=1)
	for(j in unique(temp$Year)){
		Z <- subset(temp, temp$Year == j)
		# create empty plotting window
		plot(1, type="n", xlab="", ylab="", 
			xlim=c(0.9*min(temp[,1]), 1.1*max(temp[,1])), 
			ylim=c(0, 0.1))
		# add text
		text(1.15*min(temp[,1]),0.093,unique(Z$Year), font=2, cex=1.5)

		#store CU names for next for loop
		CUs = unique(Z$CU)
		## add CDF lines
		for(i in 1:length(CUs)){
			#subset dataframe
			temp.data <- subset(Z, Z$CU == CUs[i])
			#obtain kernal density estimates
			dens = density(temp.data[,1], na.rm=T, bw="nrd")
			#plot dens
			lines(dens, lty=1, col=palette$col[palette$CU==CUs[i]], lwd=2)
		}

		## add legend
		if (Z$Year == "2014") {
			legend("topright", legend=palette$CU, lwd=2, col=palette$col)
		}	

		## add axis label
		if (Z$Year == "2016") {
		mtext(side=1, text=expression(paste("Duration of Migration in Strait")), 
		cex=1, las=1, outer=TRUE, line=1)
		mtext(side=2, text=expression(paste("Density")), 
		cex=1, las=3, outer=TRUE, line=1)
		}
	}

dev.off()

### Use this legend code inside i loop if you want a year-specific list of CUs
# legend("topright", legend=palette$CU[palette$CU %in% CUs], lwd=2,
# 			col=palette$col[palette$CU %in% CUs])

pdf("EntryDateKernel.pdf", width=6.65, height=6.65)

temp2 <- JSdat[,c("EntryDate","CU","Year")]

par(mfrow=c(3,1), mar=c(1.75,1.25,1.5,1)+0.1, oma=c(3,3,0,0))
	for(j in unique(temp2$Year)){
		Z <- subset(temp2, temp$Year == j)
		# create empty plotting window
		plot(1, type="n", xlab="", ylab="", 
			xlim=c(0.9*min(temp2[,1]), 1.1*max(temp2[,1])), ylim=c(0, 0.14),
			axes=TRUE)

		# add text
		text(0.92*min(temp2[,1]),0.13,unique(Z$Year), font=2, cex=1.5)

		#store CU names for next for loop
		CUs = unique(Z$CU)
		## add CDF lines
		for(i in 1:length(CUs)){
			#subset dataframe
			temp.data <- subset(Z, Z$CU == CUs[i])
			#obtain kernal density estimates
			dens = density(temp.data[,1], na.rm=T, bw="nrd")
			#plot dens
			lines(dens, lty=1, col=palette$col[palette$CU==CUs[i]], lwd=2)
		}

		## add legend
		if (Z$Year == "2014") {
			legend("topright", legend=palette$CU, lwd=2, col=palette$col)
		}	

		## add axis label
		if (Z$Year == "2016") {
		mtext(side=1, text=expression(paste("Ocean Entry Date")), 
		cex=1, las=1, outer=TRUE, line=1)
		mtext(side=2, text=expression(paste("Density")), 
		cex=1, las=3, outer=TRUE, line=1)
		}
}

dev.off()


########################################################
########################################################
########################################################


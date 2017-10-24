##################################################################
##################################################################
######## Johnstone Strait - Kernel Density Estimators 
###### Examine distribution of midday capture timing and duration
###### of residence for Fraser River sockeye salmon 
###### June 2017
##################################################################
##################################################################


rm(list=ls())



setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")


######## LIBRARIES ########

require(HDMD) #pairwise.mahalanobis function
require(cluster) #agnes function
require(vegan)



JSdat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))



####### DATA CLEAN #######

## Remove extra columns
d <- JSdat[,c("FishNumber","Year","CU","Age","TotalCount","EntryDate",
	"EntryFL","BodyGrowth","MigRate","MidDay")]


d$CU <- as.factor(d$CU)
d$Year <- as.factor(d$Year)


## Make palette dataframe of colors for plotting based on region
## Gradient colors
# palette <- data.frame(CU=unique(d$CU), col=c("#d73027","#f46d43","#fdae61",
# 	"#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850"))
## Discrete colors
palette <- data.frame(CU=unique(d$CU), col=c("#a6cee3","#1f78b4","#b2df8a",
	"#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"))
palette$col<-as.character(palette$col)


####### WRITE FUNCTION #######

setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Figures")


### Loop to subset dataframe based on factor levels and 
## plot multiple density plots
pdf("CountKernel.pdf", width=6.65, height=6.65)

temp <- d[,c("TotalCount","CU","Year")]

par(mfrow=c(3,1), mar=c(1.75,1.25,1.5,1)+0.1, oma=c(3,3,0,0))
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
		mtext(side=1, text=expression(paste("Days in Strait")), 
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

temp2 <- d[,c("EntryDate","CU","Year")]

par(mfrow=c(3,1), mar=c(1.75,1.25,1.5,1)+0.1, oma=c(3,3,0,0))
	for(j in unique(temp2$Year)){
		Z <- subset(temp2, temp$Year == j)
		# create empty plotting window
		plot(1, type="n", xlab="", ylab="", 
			xlim=c(0.9*min(temp2[,1]), 1.1*max(temp2[,1])), 
			ylim=c(0, 0.15))
		# add text
		text(0.92*min(temp2[,1]),0.14,unique(Z$Year), font=2, cex=1.5)

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



############## MULTIVARIATE ANALYSES ##############
#### For each year, calculate cluster analysis, NMDS then CAP (approximately
### equivalent to RDA) to look at multivariate dispersion of migratory traits among CUs


####### Aggolmerative Cluster Analysis ####### 

### Use pairwise Mahalanobis which allows groups to be easily recognized
## and accounts for original correlation structure of the variables
year = unique(JSdat$Year)


pdf("MigTraitCluster.pdf", width=6.65, height=6.65)

par(mfrow=c(2,2), mar=c(0,1.25,1.5,1)+0.1, oma=c(1,3,0,0))

for(i in 1:length(year)){
	temp <- subset(JSdat, JSdat$Year == year[i])
	
	cu = matrix(temp$CU)
	cu = t(cu[,1])

	mig_traits = c("TotalCount","BodyGrowth","EntryFL","MidDay") 
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
	plot(cluster,which.plots=2,main="") #plot dendrogram
	title(year[i])

	## store and calculate PERMANOVAs (note these use Euclidean
	## distances, not Mahalanobis) 
	m1 <- adonis(mig_traits ~ as.factor(cu), permutations=999, method="euclidean")
	assign(paste("mod", year[i], sep="."), m1)
	print(m1) 
}

mtext(side=2, text=expression(paste("Mahalanobis Distance")), las=3, outer=TRUE,
 line=1)

dev.off()

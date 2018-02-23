##################################################################
##################################################################
######## Johnstone Strait - Data Explore 
###### Examine raw data for DIPS survey 2014-2016
###### June 2017
##################################################################
##################################################################



rm(list = ls(all=TRUE)); #Remove all the objects in the memory


setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")



require(plyr)
require(fossil)
require(ggplot2)
require(vegan)
require(MASS)
require(scales) 
require(reshape2)
require(moments)
require(car)
require(HDMD) #pairwise.mahalanobis function
require(cluster) #agnes function


JSdat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))

#######################################################
######### DATA CLEAN #########

JSdat <- JSdat[,-1]


### Convert to factors
JSdat$Year <- as.factor(JSdat$Year)
JSdat$Stock <- as.factor(JSdat$Stock)
JSdat$CU <- as.factor(JSdat$CU)


## Exclude CUs with less than 10 individuals
# keep <- levels(JSdat$CU)[table(JSdat$CU)>8]
# tempJS <- JSdat[JSdat$CU %in% keep,]




#######################################################
######### DATA EXPLORE #########
table(JSdat$Stock, JSdat$Year)
table(JSdat$CU, JSdat$Year)
ddply(JSdat, .(as.factor(Year)), summarize, startSampling=min(JulianDate), endSampling=max(JulianDate),
	startEntry=min(EntryDate), endEntry=max(EntryDate), startCount=min(TotalCount), endCount=max(TotalCount))




########### SHIP DATA EXPLORE ###########



### Year specific rates
pdf("PooledBoxPlots.pdf", height=4, width=8)

par(mfrow=c(2,4), mar=c(3,3,0,0)+0.1, mgp=c(2,0.5,0))
plot(JulianDate ~ Year, JSdat)
plot(ShipFL ~ Year, JSdat)
plot(TotalCount ~ Year, data=JSdat)
plot(EntryDate ~ Year, data=JSdat)
plot(MarGrowth ~ Year, data=JSdat)
plot(Week1 ~ Year, data=JSdat)
plot(WeekF ~ Year, data=JSdat)
plot(BodyGrowth ~ Year, data=JSdat)

dev.off()


### Year and stock specific box plots 
pdf("StockSpecificPlots.pdf", height=4, width=8)

ggplot(JSdat, aes(x=Year, y=TotalCount)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=EntryFL)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=MidDay)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=BodyGrowth)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=MigRate)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=MigRateBLS)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=EntryDate)) +
	geom_boxplot() +
	facet_grid(.~Stock)

ggplot(JSdat, aes(x=Year, y=TotalCount)) +
	geom_boxplot() +
	facet_grid(.~Stock)

dev.off()


### Year and stock specific box plots 
pdf("CUSpecificPlots.pdf", height=4, width=8)

ggplot(JSdat, aes(x=Year, y=TotalCount)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=EntryFL)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=MidDay)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=BodyGrowth)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=MigRate)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=MigRateBLS)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=EntryDate)) +
	geom_boxplot() +
	facet_grid(.~CU)

ggplot(JSdat, aes(x=Year, y=TotalCount)) +
	geom_boxplot() +
	facet_grid(.~CU)

dev.off()

# JSdat[JSdat$MigRate > 50,]

##### Cumulative catch frequencies by stock
pdf("CatchFrequencies.pdf", height=4, width=8)

ggplot(JSdat, aes(JulianDate, colour = CU)) + 
	stat_ecdf() +
	facet_grid(.~Year)

ggplot(JSdat, aes(EntryDate, colour = CU)) + 
	stat_ecdf() +
	facet_grid(.~Year)

dev.off()




### Determine whether significant differences in traits among pops to justify
## using CUs with simple linear models
ED.stock <- lm(EntryDate ~ Year+Stock, data=JSdat)
ED.stock2 <- lm(EntryDate ~ Year:Stock, data=JSdat)
ED.cu <- lm(EntryDate ~ Year+CU, data=JSdat)
ED.cu2 <- lm(EntryDate ~ Year:CU, data=JSdat)
AIC(ED.stock, ED.stock2, ED.cu, ED.cu2)

ES.stock <- lm(EntryFL ~ Year+Stock, data=JSdat)
ES.stock2 <- lm(EntryFL ~ Year:Stock, data=JSdat)
ES.cu <- lm(EntryFL ~ Year+CU, data=JSdat)
ES.cu2 <- lm(EntryFL ~ Year:CU, data=JSdat)
AIC(ES.stock, ES.stock2, ES.cu, ES.cu2)

MR.stock <- lm(MigRate ~ Year+Stock, data=JSdat)
MR.stock2 <- lm(MigRate ~ Year:Stock, data=JSdat)
MR.cu <- lm(MigRate ~ Year+CU, data=JSdat)
MR.cu2 <- lm(MigRate ~ Year:CU, data=JSdat)
AIC(MR.stock, MR.stock2, MR.cu, MR.cu2)

### CUs supported for all 



### Phenology histograms
ggplot(JSdat, aes(x=MidDay, fill=CU)) +
      geom_histogram(binwidth=3, alpha=.5, position="identity") + 
      # geom_vline(data=alldat, aes(xintercept=mean(length), colour="Overall Mean"),
      #       colour="black", linetype="dashed", size=1) + 
      facet_wrap(~Year, scales="free")

ggplot(JSdat, aes(x=TotalCount, fill=CU)) +
      geom_histogram(binwidth=3, alpha=.5, position="identity") + 
      # geom_vline(data=alldat, aes(xintercept=mean(length), colour="Overall Mean"),
      #       colour="black", linetype="dashed", size=1) + 
      facet_wrap(~Year, scales="free")


### Kernel smoothed density plots


### Growth correlations
par(mfrow=c(1,1))
plot(BodyGrowth ~ MigRate, data=JSdat, pch=16, col=as.numeric(factor(Year)))
plot(OtoGrowth ~ MigRate, data=JSdat, pch=16, col=as.numeric(factor(Year)))
plot(BodyGrowth ~ MigRateBLS, data=JSdat, pch=16, col=as.numeric(factor(Year)))
plot(OtoGrowth ~ MigRateBLS, data=JSdat, pch=16, col=as.numeric(factor(Year)))
plot(BodyGrowth ~ EntryDate, data=JSdat)
plot(OtoGrowth ~ EntryDate, data=JSdat)
plot(OtoGrowth ~ EntryFL, data=JSdat)
plot(WeekF ~ EntryFL, data=JSdat)


### Residency correlations
par(mfrow=c(1,1))
plot(TotalCount ~ EntryDate, data=JSdat)
plot(TotalCount ~ EntryFL, data=JSdat)
plot(TotalCount ~ BodyGrowth, data=JSdat)

par(mfrow=c(1,1))
plot(MidDay ~ EntryDate, data=JSdat)
plot(MidDay ~ EntryFL, data=JSdat)
plot(MidDay ~ BodyGrowth, data=JSdat)


### Final week 
plot(WeekF ~ EntryDate, data=JSdat)
plot(Week1 ~ EntryDate, data=JSdat)


### Mig rate correlations 
plot(log(MigRate) ~ EntryFL, data=JSdat)
plot(log(MigRate) ~ EntryDate, data=JSdat)
plot(log(MigRate) ~ BodyGrowth, data=JSdat)



#### Summary entry size
entrysizeTable <- ddply(JSdat, .(CU, Year), summarize, n=length(EntryFL),
	CV=sd(EntryFL)/mean(EntryFL), mean=mean(EntryFL), skew=skewness(EntryFL),
	normal=round(shapiro.test(EntryFL)$p.value, 3))
entrysizeTable


ddply(JSdat[JSdat$Age==1,], .(Year), summarize, n=length(EntryFL),
	CV=sd(EntryFL)/mean(EntryFL), mean=mean(EntryFL), skew=skewness(EntryFL),
	normal=round(shapiro.test(EntryFL)$p.value, 3))

#### Summary entry date
entrydateTable <- ddply(JSdat, .(CU, Year), summarize, n=length(EntryDate),
	CV=sd(EntryDate)/mean(EntryDate), mean=mean(EntryDate), skew=skewness(EntryDate),
	normal=round(shapiro.test(EntryDate)$p.value, 3))
entrydateTable

#### Summary ocean res
resTable <- ddply(JSdat, .(CU, Year), summarize, n=length(TotalCount),
	CV=sd(TotalCount)/mean(TotalCount), mean=mean(TotalCount), skew=skewness(TotalCount),
	normal=round(shapiro.test(TotalCount)$p.value, 3))
resTable

#### Summary ocean growth
growthTable <- ddply(JSdat, .(CU, Year), summarize, n=length(BodyGrowth),
	CV=sd(BodyGrowth)/mean(BodyGrowth), mean=mean(BodyGrowth), skew=skewness(BodyGrowth),
	normal=round(shapiro.test(BodyGrowth)$p.value, 3))
growthTable

#### Summary migrateory rate BLS
migTable <- ddply(JSdat, .(CU, Year), summarize, n=length(MigRateBLS),
	CV=sd(MigRateBLS)/mean(MigRateBLS), mean=mean(MigRateBLS), skew=skewness(MigRateBLS),
	normal=round(shapiro.test(MigRateBLS)$p.value, 3))
migTable

#### Summary migrateory rate 
migTable <- ddply(JSdat, .(CU, Year), summarize, n=length(MigRate),
	CV=sd(MigRate)/mean(MigRate), mean=mean(MigRate), skew=skewness(MigRate),
	normal=round(shapiro.test(MigRate)$p.value, 3))
migTable2

#### Are data normal? - generally a mix


#### Is variance different across CUs
LeveneTable <- ddply(tempJS, .(Year), summarize, EntryFL=leveneTest(EntryFL, CU)[1,3],
	EntryDate=leveneTest(EntryDate, CU)[1,3], Res=leveneTest(TotalCount, CU)[1,3],
	MigRate=leveneTest(MigRateBLS, CU)[1,3])
LeveneTable
## Typically yes





#######################################################################################



############## MULTIVARIATE ANALYSES ##############
#### For each year, calculate cluster analysis, NMDS then CAP (approximately
### equivalent to RDA) to look at multivariate dispersion of migratory traits among CUs


####### Aggolmerative Cluster Analysis ####### 

### Use pairwise Mahalanobis which allows groups to be easily recognized
## and accounts for original correlation structure of the variables
year = unique(JSdat$Year)


par(mfrow=c(2,2))
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
	plot(cluster,which.plots=2,main="") #plot dendrogram
	title(year[i])

	## store and calculate PERMANOVAs (note these use Euclidean
	## distances, not Mahalanobis) 
	m1 <- adonis(mig_traits ~ as.factor(cu), permutations=999, method="euclidean")
	assign(paste("mod", year[i], sep="."), m1)
	print(m1) 
}



####### NMDS ####### 

library(ecodist)

year = unique(JSdat$Year)


pal <- data.frame(CU=unique(JSdat$CU), col=c("#a1d99b","#deebf7","#fc9272",
	"#9ecae1","#de2d26","#3182bd","#31a354","#bcbddc","#756bb1"))
pal$col<-as.character(pal$col)


pdf("MigTrait_NMDS.pdf", height=7, width=7)

par(mfrow=c(2,2), mar=c(1,1,1,0.1), oma=c(1.5,1.5,1,1))
for(i in 1:length(year)){
	temp <- subset(JSdat, JSdat$Year == year[i])
	
	cu = matrix(temp$CU)
	cu = t(cu[,1])

	mig_traits_col = c("TotalCount","BodyGrowth","EntryFL",
		"EntryDate") 
	mig_traits = as.matrix(temp[,mig_traits_col])

	dm = distance(mig_traits, "mahal")

	nmds1 <- metaMDS(dm, k=2, trymax=500)

	cu <- temp$CU

	plot(nmds1$points, pch=1, main=year[i], cex=0.8)
	vectors = vf(nmds1$points, temp[,mig_traits_col], nperm=10)
	plot(vectors, len=0.1, col="red")
	
	CUs = sort(unique(temp$CU))
	
	for (i in 1:length(CUs)){
 		ordiellipse(nmds1$points, temp$CU,
 		conf=0.80, 
 		col=pal$col[pal$CU==CUs[i]],
 		show.groups=CUs[i],lwd=2)
 	}
}
### add legend to last panel
plot(nmds1$points, pch=19, col="white", axes=FALSE, xlab="", ylab="")
legend("topright", legend=pal$CU, lwd=2,
			col=pal$col)

dev.off()

####### CAP ########


mig_traits = c("TotalCount","BodyGrowth","EntryFL","EntryDate", "MidDay") 
	### NOTE: cannot include midday and entry date and calculate mahalanobis
	### distance; too strongly correlated; AC higher for entrydate so 
	### currently using those
	mig_traits = as.matrix(JSdat[,mig_traits])


### Fit constrained ordination after "removing" year effect
JScap <- capscale(mig_traits ~ CU + Condition(Year), JSdat, dist="euclidean")
JScap
coef(JScap)
RsquareAdj(JScap)$adj.r.squared  ##0.15


### Plot and test for significant differences
plot(JScap, scaling=2, display=c("sp","lc","cn"))
js.sc <- scores(JScap, choices=1:2, scaling=2, display="sp")
arrows(0, 0, js.sc[,1], js.sc[,2], length=0, lty=1, col="red")

# anova(JScap)

#############################
#############################
#############################


############################################
############################################


####### ANALYSIS OF MULTIVARIATE DISPERSIONS ####### 
JSmat <- JSdat[,c("TotalCount","BodyGrowth","EntryFL","EntryDate","MidDay")]
JSmat <- decostand(JSmat, "standardize")

mig.dist <- vegdist(JSmat, method="euclidean")
mod <- betadisper(mig.dist, JSdat$CU)


plot(TukeyHSD(mod))
permutest(mod, pairwise=TRUE)


############################################
############################################
############################################
############################################
############################################
############################################



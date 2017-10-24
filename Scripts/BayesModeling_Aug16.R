##################################################################
##################################################################
######## Johnstone Strait - Hierarchical Modeling 
###### Use top ranked models to explore individual, population and 
###### year scale effects on phenology timing and duration; different
###### from no date script because of random slopes for age; different
###### from July script because of added age-2 effect
###### August 2017
##################################################################
##################################################################


rm(list=ls())


setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")



require(rethinking)
require(car)
require(gplots)

JSdat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))



####### DATA CLEAN #######

JSdat$Age <- as.factor(JSdat$Age)


#### Remove low sample size CUs and re-run models
## Trim at 10
# temp14 <- JSdat[JSdat$Year == "2014",]
# temp15 <- JSdat[JSdat$Year == "2015" & !JSdat$CU == "Francois-Fraser-S" 
# 	& !JSdat$CU == "Lillooet-Harrison-L" & !JSdat$CU == "Nadina-Francois-ES",]
# temp16 <- JSdat[JSdat$Year == "2016",]

## Trim at 20
# temp14 <- JSdat[JSdat$Year == "2014" & !JSdat$CU == "Pitt-ES",]
# temp15 <- JSdat[JSdat$Year == "2015" & !JSdat$CU == "Francois-Fraser-S" 
# 	& !JSdat$CU == "Lillooet-Harrison-L" & !JSdat$CU == "Nadina-Francois-ES",]
# temp16 <- JSdat[JSdat$Year == "2016" & !JSdat$CU == "Francois-Fraser-S" 
# 	& !JSdat$CU == "Lillooet-Harrison-L" & !JSdat$CU == "Nadina-Francois-ES",]

# JSdat <- rbind(temp14, temp15, temp16)



## Make index variables
JSdat$cu_id <- coerce_index(JSdat$CU)
JSdat$yr_id <- coerce_index(JSdat$Year)
JSdat$age_id <- 0
JSdat$age_id[JSdat$Age=="0"] <- 1
JSdat$age2_id <- 0
JSdat$age2_id[JSdat$Age=="2"] <- 1

## Transform explanatory variables
JSdat$fl_t <- (JSdat$EntryFL - mean(JSdat$EntryFL))/sd(JSdat$EntryFL) 
JSdat$date_t <- (JSdat$EntryDate - mean(JSdat$EntryDate))/sd(JSdat$EntryDate)
JSdat$growth_t <- (JSdat$BodyGrowth - mean(JSdat$BodyGrowth))/sd(JSdat$BodyGrowth)
JSdat$river_t <- (JSdat$RiverDist - mean(JSdat$RiverDist))/sd(JSdat$RiverDist)


## Remove extra columns
d <- JSdat[,c("FishNumber","yr_id","cu_id","age_id","age2_id","TotalCount","EntryDate",
	"fl_t","date_t","growth_t","river_t")]


# plot(TotalCount ~ fl_t, data=d)


### Average differences across age restricted to stocks with age-2
temp <- JSdat[JSdat$Stock == "Chilko" | JSdat$Stock == "Birkenhead",]
ddply(temp, .(as.factor(Age)), summarize, meanCount=mean(TotalCount), sdCount=sd(TotalCount), meanBodyGrowth=mean(BodyGrowth), 
		sdBodyGrowth=sd(BodyGrowth), meanEntryFL=mean(EntryFL), sdEntryFL=sd(EntryFL), meanEntryDate=mean(EntryDate), sdEntryDate=var(EntryDate))

################################################################
################################################################
################################################################


### Total count models
m1.dgs1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bs_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		c(bd, bg, bs) ~ dnorm(0,3),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)


### Entry date models
m3.rs1 <- map2stan(
	alist(
		#likelihood
		EntryDate ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A +  BR*river_t + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BR <- br + br_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bs_yr, br_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		c(bs,br) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)


### Check convergence
precis(m1.dgs1, depth=2)
precis(m3.rs1, depth=2)


plot(precis(m1.dgs1, depth=2, pars=c("bd","bg","bs","bg_cu","bs_cu","bd_cu")))
plot(precis(m1.dgs1, depth=2, pars=c("bd","bg","bs","bg_yr","bs_yr","bd_yr")))
################################################################
################################################################
################################################################



############## LINEAR MODEL FIGURES ##############



########## Figure 5 - estimates of age-1 and age-2 effect sizes (two panel)


########## Calculate age effects 

##### Create dataframe with average values and HPDI for entry date age effects
DateAge <- data.frame(mean=coef(m3.rs1)[c(53,54)], l=rep(0,2), h=rep(0,2))
## Matrix of year specific posterior samples from effects distributions
DateMeanAge <- matrix(as.numeric(
	c(post2$ba,post2$bb)), ncol=2)

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(DateAge)){
	low=HPDI(DateMeanAge[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(DateMeanAge[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

DateAge$l <- as.numeric(low.vec[1:2])
DateAge$h <- as.numeric(high.vec[1:2])
DateAge$Model <- as.factor("date")



##### Create dataframe with average values and HPDI for duration sigma
CountAge <- data.frame(mean=coef(m1.dgs1)[c(77,78)], l=rep(0,2), h=rep(0,2))
## Matrix of year specific posterior samples from effects distributions
CountMeanAge <- matrix(as.numeric(
	c(post1$ba,post1$bb)), ncol=2)

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(CountAge)){
	low=HPDI(CountMeanAge[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(CountMeanAge[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

CountAge$l <- as.numeric(low.vec[1:2])
CountAge$h <- as.numeric(high.vec[1:2])
CountAge$Model <- as.factor("count")


ageDF <- rbind(DateAge,CountAge)
ageDF <- ageDF[c(1,3,2,4),]


###### Figure
age_col = c("#bcbddc","#bcbddc","#756bb1","#756bb1")


png("AgeEstimate.png", width=3.5, height=3.5, units="in", res=900)

par(mar=c(2.4,3,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3)

plotCI(c(1.35,1.65,2.85,3.15), ageDF$mean, ui=(ageDF$h), li=(ageDF$l), ylim=c(-10,50), xlim=c(0.8,3.65), xlab="",
	ylab="",axes=FALSE, pch=c(24,21,24,21), pt.bg=age_col, cex=2, sfrac=0)
axis(1, tick=T, at=c(-1,1.5,3,5), lwd=2, lwd.ticks=1, cex.axis=1, 
	labels=c("","Age-0","Age-2",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-15,-10,15,45,55),
  labels=c("","-10","15","45",""), las=1)
mtext(side=2, line=1.6, text=expression(paste("Age Effect (Days Difference)" )), 
	outer=FALSE, cex=1)
abline(h=0, lwd=2, lty=2)

dev.off()


######################################################################
######################################################################



########## Figure 6 - Estimates of slope parameters for entry date models
 

######## Create data frame with hyper- and yrspecific means for each variable


### 1) Calculate slope parameters for entry date (entry size and river distance
## covariate)
### Calculate slope parameters for total count models across all covariates
post2 <- extract.samples(m3.rs1)


### Create matrix of distributions
muDateBetas <- matrix(as.numeric(
	rbind(post2$bs,post2$br)), nrow=2)
rownames(muDateBetas) = c("bs","br")

### Calculate HPDIs for each beta
datePI <- apply(muDateBetas, 1, function(x) HPDI(x, prob=0.9))

### Combine
muDF <- data.frame(mean=coef(m3.rs1)[c(55,56)], l=datePI[1,], h=datePI[2,])


## compare
muDF
precis(m3.rs1) ### looks good


######### Yr specific means

##### Create dataframe with average values and HPDI for each year specific 
#### coefficient estimate
yrDateBetas <- data.frame(mean=coef(m3.rs1)[c(46:48,43:45)], l=rep(0,6), h=rep(0,6))
## Matrix of year specific posterior samples from effects distributions
yrDateMeanBetas <- matrix(as.numeric(c(post2$bs_yr,post2$br_yr)), ncol=6)
colnames(yrDateMeanBetas) = (c(rep("bs_yr",3),rep("br_yr",3)))

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(yrDateBetas)){
	low=HPDI(yrDateMeanBetas[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(yrDateMeanBetas[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

yrDateBetas$l <- as.numeric(low.vec[1:6])
yrDateBetas$h <- as.numeric(high.vec[1:6])


### Compare new data frame to precis (looks good!)
yrDateBetas
precis(m3.rs1, depth=2, pars=c("bs_yr","br_yr"))



########### Figure 

setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/Figures")


#### Common graphical parameters
ylimits = c(-2,5)

#### Entry size coefficient colors
col_bs = c("#eff3ff","#bdd7e7","#6baed6")
col_bsMean = "#2171b5"
col_br = c("#fee5d9","#fcae91","#fb6a4a")
col_brMean = "#cb181d"

#### Layout
mat<-matrix(1, 100, 100)
mat[1:100,26:100]<-2
# mat[1:50,26:100]<-2
# mat[51:100,1:25]<-3
# mat[51:100,26:100]<-4

png("EntryDateSlope.png", width=3.5, height=2, units="in", res=900)

par(mar=c(2.4,2.4,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3)

layout(mat)

######## REMOVE RIVER EFFECTS FOR NOW 
### Top row: river distance effects on entry date
# plotCI(1,muDF$mean[2], ui=(muDF$h[2]), li=(muDF$l[2]), ylim=ylimits, xlab="",
# 	ylab="",axes=FALSE, pch=24, pt.bg=col_brMean, cex=2, sfrac=0)
# axis(1, tick=T, at=c(0,1,3), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","Hypermean",""), tcl=0) 
# axis(2, tick=T, lwd=2, cex.axis=1, at=c(-4,0,4,7),
#   labels=c("","0","4",""), las=1)
# mtext(side=2, line=1.2, text=expression(paste("River Distance Effect" )), 
# 	outer=FALSE, cex=0.8)
# abline(h=0, lty=2)
# text(1, 4.6, "a)", font=1.5, cex=1.5) 

# deviations from mean effect
# plotCI(c(1:3), yrDateBetas$mean[4:6], ui=c(yrDateBetas$h[4:6]), li=c(yrDateBetas$l[4:6]), 
# 	xlim=c(0.8,3.2), ylim=ylimits, xlab="", ylab="", axes=FALSE, pch=24, pt.bg=col_br, cex=2, sfrac=0)
# axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","2014","2015","2016",""), 
# 	tcl=0)
# abline(h=0, lty=2)

### Bottom row: size effects on entry date
plotCI(1,muDF$mean[1], ui=(muDF$h[1]), li=(muDF$l[1]), ylim=ylimits, xlab="",
	ylab="",axes=FALSE, pch=24, pt.bg=col_bsMean, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,3), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","Hypermean",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-4,0,4,7),
  labels=c("","0","4",""), las=1)
mtext(side=2, line=1.2, text=expression(paste("Entry Size Effect" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
# text(1, 4.6, "b)", font=1.5, cex=1.5) 

# deviations from mean effect
plotCI(c(1:3), yrDateBetas$mean[1:3], ui=c(yrDateBetas$h[1:3]), li=c(yrDateBetas$l[1:3]), 
	xlim=c(0.8,3.2), ylim=ylimits, xlab="", ylab="", axes=FALSE, pch=24, pt.bg=col_bs, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","2014","2015","2016",""), 
	tcl=0)
abline(h=0, lty=2)


dev.off()


######################################################################
######################################################################






########## Figure 7 - four panel dot plot of effect sizes for all covariates on 
##### total counts


### Calculate slope parameters for total count models across all covariates
post1 <- extract.samples(m1.dgs1)


### Create matrix of distributions
muCountBetas <- matrix(as.numeric(
	rbind(post1$bs,post1$bd,post1$bg)), nrow=3)
rownames(muCountBetas) = c("bs","bd","bg")

### Calculate HPDIs for each beta
countPI <- apply(muCountBetas, 1, function(x) HPDI(x, prob=0.9))

### Combine
muDF <- data.frame(mean=coef(m1.dgs1)[c(76,74,75)], l=countPI[1,], h=countPI[2,])


## compare
muDF
precis(m1.dgs1)


######### Yr specific means

##### Create dataframe with average values and HPDI for each year specific 
#### coefficient estimate
yrCountBetas <- data.frame(mean=coef(m1.dgs1)[c(64:69,61:63)], l=rep(0,9), h=rep(0,9))
## Matrix of year specific posterior samples from effects distributions
yrCountMeanBetas <- matrix(as.numeric(
	c(post1$bs_yr,post1$bd_yr,post1$bg_yr)), ncol=9)
colnames(yrCountMeanBetas) = (c(rep("bs_yr",3),rep("bd_yr",3),rep("bg_yr",3)))

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(yrCountBetas)){
	low=HPDI(yrCountMeanBetas[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(yrCountMeanBetas[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

yrCountBetas$l <- as.numeric(low.vec[1:9])
yrCountBetas$h <- as.numeric(high.vec[1:9])


### Compare new data frame to precis (looks good!)
yrCountBetas
precis(m1.dgs1, depth=2, pars=c("bs_yr","bd_yr","bg_yr"))



########### Figure 

#### Common graphical parameters
ylimits = c(-8,10)

#### Entry size coefficient colors
col_bs = c("#eff3ff","#bdd7e7","#6baed6")
col_bd = c("#edf8e9","#bae4b3","#74c476")
col_bg = c("#feedde","#fdbe85","#fd8d3c")
col_Means = c("#2171b5","#238b45","#d94701")

#### Layout
mat<-matrix(1, 100, 100)
mat[1:50,51:100]<-2
mat[51:100,1:50]<-3
mat[51:100,51:100]<-4


png("CountSlopes.png", width=5, height=5, units="in", res=900)

par(mar=c(2.4,3,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3)

layout(mat)

##overall means plot
plotCI(c(1:3),muDF$mean, ui=(muDF$h), li=(muDF$l), xlim=c(0.8,3.2), ylim=ylimits, xlab="",
	ylab="",axes=FALSE, pch=21, pt.bg=col_Means, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, 
	labels=c("","Entry Size","Entry Date","Marine Growth",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-10,-5,0,5,10,12),
  labels=c("","-5","0","5","10",""), las=1)
mtext(side=2, line=1.6, text=expression(paste("Mean Effect Sizes" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
text(0.9, 9, "a)", font=1.5, cex=1.5) 

##entry size
plotCI(c(1:3), yrCountBetas$mean[1:3], ui=(yrCountBetas$h[1:3]), li=(yrCountBetas$l[1:3]), 
	xlim=c(0.8,3.2), ylim=ylimits, xlab="", ylab="",axes=FALSE, pch=21, pt.bg=col_bs, 
	cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, 
	labels=c("","2014","2015","2016",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-10,-5,0,5,10,12),
  labels=c("","-5","0","5","10",""), las=1)
mtext(side=2, line=1.6, text=expression(paste("Year-Specific Size Effect" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
text(0.9, 9, "b)", font=1.5, cex=1.5)

##entry date
plotCI(c(1:3), yrCountBetas$mean[4:6], ui=(yrCountBetas$h[4:6]), li=(yrCountBetas$l[4:6]), 
	xlim=c(0.8,3.2), ylim=ylimits, xlab="", ylab="",axes=FALSE, pch=21, pt.bg=col_bd, 
	cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, 
	labels=c("","2014","2015","2016",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-10,-5,0,5,10,12),
  labels=c("","-5","0","5","10",""), las=1)
mtext(side=2, line=1.6, text=expression(paste("Year-Specific Date Effect" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
text(0.9, 9, "c)", font=1.5, cex=1.5) 

##marine growth
plotCI(c(1:3), yrCountBetas$mean[7:9], ui=(yrCountBetas$h[7:9]), li=(yrCountBetas$l[7:9]), 
	xlim=c(0.8,3.2), ylim=ylimits, xlab="", ylab="",axes=FALSE, pch=21, pt.bg=col_bg, 
	cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, 
	labels=c("","2014","2015","2016",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-10,-5,0,5,10,12),
  labels=c("","-5","0","5","10",""), las=1)
mtext(side=2, line=1.6, text=expression(paste("Year-Specific Growth Effect" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
text(0.9, 9, "d)", font=1.5, cex=1.5) 

 
dev.off()


par(mfrow=c(2,2))
plot(TotalCount ~ fl_t, data=d[d$yr_id=="2",], xlim=c(-2,2))
plot(TotalCount ~ fl_t, data=d[d$yr_id=="3",])

plot(fl_t~yr_id, data=d)
######################################################################
######################################################################




######## PRESENTATION FIGURE - SIZE EFFECTS ONLY 


png("CountSizeSlope.png", width=3.5, height=2, units="in", res=900)

par(mar=c(2.4,2.4,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3)

mat<-matrix(1, 100, 100)
mat[1:100,26:100]<-2

layout(mat)

plotCI(1,muDF$mean[1], ui=(muDF$h[1]), li=(muDF$l[1]), ylim=c(-6,7), xlab="",
	ylab="",axes=FALSE, pch=21, pt.bg=col_bsMean, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,3), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","Hypermean",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-8,-4,0,4,8),
  labels=c("","-4","0","4",""), las=1)
mtext(side=2, line=1.2, text=expression(paste("Entry Size Effect" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
# text(1, 4.6, "b)", font=1.5, cex=1.5) 

# deviations from mean effect
plotCI(c(1:3), yrCountBetas$mean[1:3], ui=c(yrCountBetas$h[1:3]), li=c(yrCountBetas$l[1:3]), 
	xlim=c(0.8,3.2), ylim=c(-6,7), xlab="", ylab="", axes=FALSE, pch=21, pt.bg=col_bs, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","2014","2015","2016",""), 
	tcl=0)
abline(h=0, lty=2)


dev.off()




png("DateDistSlope.png", width=3.5, height=2, units="in", res=900)

par(mar=c(2.4,2.4,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3)

layout(mat)

## Top row: river distance effects on entry date
plotCI(1,muDF$mean[2], ui=(muDF$h[2]), li=(muDF$l[2]), ylim=ylimits, xlab="",
	ylab="",axes=FALSE, pch=24, pt.bg=col_brMean, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,3), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","Hypermean",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-4,0,4,7),
  labels=c("","0","4",""), las=1)
mtext(side=2, line=1.2, text=expression(paste("River Distance Effect" )), 
	outer=FALSE, cex=0.8)
abline(h=0, lty=2)
text(1, 4.6, "a)", font=1.5, cex=1.5) 

## deviations from mean effect
plotCI(c(1:3), yrDateBetas$mean[4:6], ui=c(yrDateBetas$h[4:6]), li=c(yrDateBetas$l[4:6]), 
	xlim=c(0.8,3.2), ylim=ylimits, xlab="", ylab="", axes=FALSE, pch=24, pt.bg=col_br, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, labels=c("","2014","2015","2016",""), 
	tcl=0)
abline(h=0, lty=2)

dev.off()





######################################################################
######################################################################





########## Figure 8 - dot plot of distribution of average and CU-specific mean values
######## (i.e. intercepts)



########### Calculate parameters for entry date

post2 <- extract.samples(m3.rs1)

### Extract average mean
meanDate <- coef(m3.rs1)[52]
muDate <- post2$a
muDatePI <- HPDI(muDate, prob=0.9)


##### Create dataframe with average values and HPDI for each year
### Extract vector of year averages
yrDates <- data.frame(mean=coef(m3.rs1)[49:51], l=rep(0,3), h=rep(0,3))
yrMeans <- post2$a_yr

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(yrDates)){
	low=HPDI(yrMeans[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(yrMeans[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

yrDates$l <- as.numeric(low.vec[1:3])
yrDates$h <- as.numeric(high.vec[1:3])

##### Create dataframe with average values and HPDI for each CU
### Extract vector of CU averages
cuDates <- data.frame(mean=coef(m3.rs1)[28:36], l=rep(0,9), h=rep(0,9))
cuMeans <- post2$a_cu

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(cuDates)){
	low=HPDI(cuMeans[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(cuMeans[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

cuDates$l <- as.numeric(low.vec[1:9])
cuDates$h <- as.numeric(high.vec[1:9])

#### Check estimates
yrDates
cuDates
precis(m3.rs1, depth=2, pars=c("a","a_yr","a_cu"))

#######################################


########### Calculate parameters for days resident

post1 <- extract.samples(m1.dgs1)

### Extract average mean
meanCount <- coef(m1.dgs1)[73]
muCount <- post1$a
muPI <- HPDI(muCount, prob=0.9)


##### Create dataframe with average values and HPDI for each year
### Extract vector of year averages
yrCounts <- data.frame(mean=coef(m1.dgs1)[70:72], l=rep(0,3), h=rep(0,3))
yrMeans <- post1$a_yr

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(yrCounts)){
	low=HPDI(yrMeans[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(yrMeans[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

yrCounts$l <- as.numeric(low.vec[1:3])
yrCounts$h <- as.numeric(high.vec[1:3])

##### Create dataframe with average values and HPDI for each CU
### Extract vector of CU averages
cuCounts <- data.frame(mean=coef(m1.dgs1)[46:54], l=rep(0,9), h=rep(0,9))
cuMeans <- post1$a_cu

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(cuCounts)){
	low=HPDI(cuMeans[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(cuMeans[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

cuCounts$l <- as.numeric(low.vec[1:9])
cuCounts$h <- as.numeric(high.vec[1:9])


### Compare new data frame to precis (looks good!)
cuCounts
precis(m1.dgs1, depth=2, pars=c("a","a_cu"))



###### FIGURE  ######
#### Construct six panel matrix to plot overall average, year specific, then CU-specific means 


## Make palette dataframe of colors for plotting based on region
## Discrete colors
temp <- unique(JSdat[c("CU", "cu_id")])
temp <- temp[order(temp$cu_id),] # reorder based on CU id
palette <- data.frame(CU=temp$CU, cu_id=temp$cu_id, col=c("#a6cee3","#1f78b4","#b2df8a",
	"#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"))
palette$col<-as.character(palette$col)


### Plotting area
mat<-matrix(1, 100, 100)
mat[1:50,16:35]<-2
mat[1:50,36:100]<-3
mat[51:100,1:16]<-4
mat[51:100,16:35]<-5
mat[51:100,36:100]<-6


#### Common graphical parameters
ylimits = c(-12,12)
year_col = c("white", "lightgrey", "darkgrey")



png("InterceptEstimates.png", width=7, height=5, units="in", res=900)

par(mar=c(2.4,3.5,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3, cex.axis=1)

layout(mat)
# mean entry date
plotCI(1,meanDate, ui=(muDatePI[2]), li=(muDatePI[1]), ylim=c(meanDate-12, meanDate+12), xlab="",
	ylab="",axes=FALSE, pch=24, pt.bg="black", cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,3), lwd=2, lwd.ticks=1, labels=c("","Hypermean",""), tcl=0) 
axis(2, tick=T, lwd=2, at=c(100,115,125,135,140),
  labels=c("","115","125","135",""), las=1)
mtext(side=2, line=2.1, text=expression(paste("Entry Date (day of year)" )), 
	outer=FALSE, cex=0.9)
text(0.8, 0.99*(meanDate+12), "a)", font=1.5, cex=1.5) 

# deviations from mean entry date
plotCI(c(1:3), yrDates$mean, ui=c(yrDates$h), li=c(yrDates$l), xlim=c(0.8,3.2), ylim=ylimits, xlab="",
	ylab="", axes=FALSE, pch=24, pt.bg=year_col, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, labels=c("","2014","2015","2016",""), 
	tcl=0) 

plotCI(c(1:9), cuDates$mean, ui=c(cuDates$h), li=c(cuDates$l), ylim=ylimits, xlab="",
	ylab="",axes=FALSE, cex=2, pch=24, pt.bg=palette$col, sfrac=0)
axis(1, tick=T, at=c(0,5,10), lwd=2, lwd.ticks=1, labels=c("","Conservation Units",""), tcl=0) 
text(1:9, -12, rep(1:9), cex=1.2)
# legend("bottomright", legend=palette$CU, pch=21, pt.bg=palette$col, ncol=3, cex=1.1, pt.cex=1.4)

# mean migration duration
plotCI(1,meanCount, ui=(muPI[2]), li=(muPI[1]), ylim=c(meanCount-12, meanCount+12), xlab="",
	ylab="",axes=FALSE, pch=21, pt.bg="black", cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,3), lwd=2, lwd.ticks=1, labels=c("","Hypermean",""), tcl=0) 
axis(2, tick=T, lwd=2, at=c(0,25,35,45,70),
  labels=c("","25","35","45",""), las=1)
mtext(side=2, line=2.1, text=expression(paste("Duration of Migration (days)" )), 
	outer=FALSE, cex=0.9) 
text(0.8, 0.98*(meanCount+12), "b)", font=1.5, cex=1.5) 

# deviations from mean entry date
plotCI(c(1:3), yrCounts$mean, ui=c(yrCounts$h), li=c(yrCounts$l), xlim=c(0.8,3.2), ylim=ylimits, xlab="",
	ylab="",axes=FALSE, pch=21, pt.bg=year_col, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, labels=c("","2014","2015","2016",""), 
	tcl=0) 

plotCI(c(1:9), cuCounts$mean, ui=c(cuCounts$h), li=c(cuCounts$l), ylim=ylimits, xlab="",
	ylab="",axes=FALSE, pch=21, pt.bg=palette$col, cex=2, sfrac=0)
axis(1, tick=T, at=c(0,5,10), lwd=2, lwd.ticks=1, labels=c("","Conservation Units",""), tcl=0) 
text(1:9, -12, rep(1:9), cex=1.2)


dev.off()


######################################################################
######################################################################



########## Figure 9 - estimates of variation (sigmas) for entry date and
######## duration models (two panel)

sig_col = c("white","white","lightgrey","lightgrey","darkgrey","darkgrey")


##### Create dataframe with average values and HPDI for entry date sigma
DateSigmas <- data.frame(mean=coef(m3.rs1)[c(57,62,58)], l=rep(0,3), h=rep(0,3))
## Matrix of year specific posterior samples from effects distributions
DateMeanSigma <- matrix(as.numeric(
	c(post2$sigma,post2$sigma_yr[,1],post2$sigma_cu[,1])), ncol=3)

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(DateSigmas)){
	low=HPDI(DateMeanSigma[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(DateMeanSigma[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

DateSigmas$l <- as.numeric(low.vec[1:3])
DateSigmas$h <- as.numeric(high.vec[1:3])
DateSigmas$Model <- as.factor("date")

precis(m3.rs1, depth=2, pars=c("sigma","sigma_cu","sigma_yr"))


##### Create dataframe with average values and HPDI for duration sigma
CountSigmas <- data.frame(mean=coef(m1.dgs1)[c(79,86,80)], l=rep(0,3), h=rep(0,3))
## Matrix of year specific posterior samples from effects distributions
CountMeanSigma <- matrix(as.numeric(
	c(post1$sigma,post1$sigma_yr[,1],post1$sigma_cu[,1])), ncol=3)

low.vec=NULL
high.vec=NULL
for(i in 1:nrow(CountSigmas)){
	low=HPDI(CountMeanSigma[,i], prob=0.9)[1] #lower bound to interval
	low.vec=c(low.vec, low)

	high=HPDI(CountMeanSigma[,i], prob=0.9)[2] #upper bound to interval
	high.vec=c(high.vec, high)
}

CountSigmas$l <- as.numeric(low.vec[1:3])
CountSigmas$h <- as.numeric(high.vec[1:3])
CountSigmas$Model <- as.factor("count")

sigmaDF <- rbind(DateSigmas,CountSigmas)
sigmaDF <- sigmaDF[c(1,4,2,5,3,6),]

###### Figure

png("SigmaEstimate.png", width=3.5, height=3.5, units="in", res=900)


par(mar=c(2.4,3,0,0)+0.1, oma=c(0,0,0,0)+0.2, mgp=c(3,0.5,0), tcl=-0.3)


plotCI(c(1.25,1.75,3.25,3.75,5.25,5.75), sigmaDF$mean, ui=(sigmaDF$h), li=(sigmaDF$l), ylim=c(0,11), xlim=c(0.8,6.2), xlab="",
	ylab="",axes=FALSE, pch=c(24,21)[as.numeric(sigmaDF$Model)], pt.bg=sig_col, cex=2, sfrac=0)
axis(1, tick=T, at=c(-1,1.5,3.5,5.5,7), lwd=2, lwd.ticks=1, cex.axis=1, 
	labels=c("","Hyperestimate","Year","CU",""), tcl=0) 
axis(2, tick=T, lwd=2, cex.axis=1, at=c(-4,0,4,8,12),
  labels=c("","0","4","8",""), las=1)
mtext(side=2, line=1.6, text=expression(paste("Sigma Estimate" )), 
	outer=FALSE, cex=1)

dev.off()


######################################################################
######################################################################














######################################################################
######################################################################

########## Figure S1  - dot plot of effect sizes of effect sizes for growth, date, size 
######## (i.e. year specific slopes)
pdf("FigS1_CUSlopes_Supp.pdf", width=6, height=6)

par(mfrow=c(2,2), mar=c(3,1,2,1), cex.main=0.9, mgp=c(2,0.5,0))
plot(precis(m3.rs1, depth=2, pars=c("bs_cu")), xlab="Effect Size", 
	main="Size Effects on Entry Date")
plot(precis(m1.dgs1, depth=2, pars=c("bs_cu")), xlab="Effect Size", 
	main="Size Effects on Migration Duration")
plot(precis(m1.dgs1, depth=2, pars=c("bd_cu")), xlab="Effect Size", 
	main="Date Effects on Migration Duration")
plot(precis(m1.dgs1, depth=2, pars=c("bg_cu")), xlab="Effect Size", 
	main="Growth Effects on Migration Duration")

dev.off()


########## Figure S2  - dot plot of effect sizes for river distance

pdf("FigS2_YrSlopesRiver_Supp.pdf", width=3, height=3)

par(mfrow=c(1,1), mar=c(3,1,3,1), cex.main=0.9, mgp=c(2,0.5,0))
plot(precis(m3.rs1, depth=2, pars=c("br_yr")), xlab="Effect Size", 
	main="River Distance Effects\non Entry Date")

dev.off()


############## CHECK PREDICTIONS ##############

########### Figure S1 - counterfactural predictions of effect sizes


### Compute counterfactual predictions for effect of size on residence
## and midday; hold growth and entry date at mean values

post <- extract.samples(m1.dgs1)

## Sample intercepts and slopes from distribution using estimates of 
## variation in cu and yr; MAKE SURE ENTERED IN CORRECT ORDER
# cu effects
a_cu_sims <- matrix(rnorm(9000, 0, post$sigma_cu[1]), 1000, 9)
bg_cu_sims <- matrix(rnorm(9000, 0, post$sigma_cu[2]), 1000, 9)
bs_cu_sims <- matrix(rnorm(9000, 0, post$sigma_cu[3]), 1000, 9)
bd_cu_sims <- matrix(rnorm(9000, 0, post$sigma_cu[4]), 1000, 9)

# yr effects
a_yr_sims <- matrix(rnorm(3000, 0, post$sigma_yr[1]), 1000, 3)
bd_yr_sims <- matrix(rnorm(3000, 0, post$sigma_yr[2]), 1000, 3)
bs_yr_sims <- matrix(rnorm(3000, 0, post$sigma_yr[4]), 1000, 3)
bg_yr_sims <- matrix(rnorm(3000, 0, post$sigma_yr[3]), 1000, 3)
ba_yr_sims <- matrix(rnorm(3000, 0, post$sigma_yr[5]), 1000, 3)


##### Predictions for date effects
date.pred <- list(
	# explanatory variable of interest
	date_t = seq(from=-4, to=4, length.out=50),
	# random effects
	cu_id = rep(1, 50),
	yr_id = rep(1, 50),
	# mean fixed effects
	fl_t = rep(mean(d$fl_t), 50),
	growth_t = rep(mean(d$growth_t), 50),
	# set age = to 1
	age_id = rep(0, 50)
	)


## Link function for intercept/slope model
link.m1.dgs1.date <- link(m1.dgs1, n=1000, data=date.pred,
	replace=list(a_cu=a_cu_sims, 
	bg_cu=bg_cu_sims, bs_cu=bs_cu_sims, bd_cu=bd_cu_sims,
	a_yr=a_yr_sims, bg_yr=bg_yr_sims, bs_yr=bs_yr_sims, bd_yr=bd_yr_sims
	))
################


##### Predictions for size effects
size.pred <- list(
	# explanatory variable of interest
	fl_t = seq(from=-2, to=5, length.out=50),
	# random effects
	cu_id = rep(1, 50),
	yr_id = rep(1, 50),
	# mean fixed effects
	date_t = rep(mean(d$date_t), 50),
	growth_t = rep(mean(d$growth_t), 50),
	# set age = to 1
	age_id = rep(0, 50)
	)


## Link function for intercept/slope model
link.m1.dgs1.size <- link(m1.dgs1, n=1000, data=size.pred,
	replace=list(a_cu=a_cu_sims, 
	bg_cu=bg_cu_sims, bs_cu=bs_cu_sims, bd_cu=bd_cu_sims,
	a_yr=a_yr_sims, bg_yr=bg_yr_sims, bs_yr=bs_yr_sims, bd_yr=bd_yr_sims
	))


## Raw data
plot(TotalCount ~ date_t, data=d, col=rangi2, pch=16, xlab="Entry Date",
	ylab="Days Resident")


## Posterior median and PI for intercept and 1 slope 
mu.median <- apply(link.m1.dgs1.date$mu, 2, median)
lines(d.pred$date_t, mu.median, lwd=2)
mu.PI <- apply(link.m1.dgs1.date$mu, 2, PI, prob=0.9)
shade(mu.PI, d.pred$date_t)


plot(TotalCount ~ fl_t, data=d, col=rangi2, pch=16, xlab="Entry Size",
	ylab="Days Resident", xlim=c(-2,5))


## Posterior median and PI for intercept and 1 slope 
mu.median <- apply(link.m1.dgs1.size$mu, 2, median)
lines(size.pred$fl_t, mu.median, lwd=2)
mu.PI <- apply(link.m1.dgs1.size$mu, 2, PI, prob=0.9)
shade(mu.PI, size.pred$fl_t)



######### Supplementary models


# m1.dgs2 <- map2stan(
# 	alist(
# 		#likelihood
# 		TotalCount ~ dnorm(mu, sigma),
		
# 		#linear models
# 		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id,
# 		A <- a + a_cu[cu_id] + a_yr[yr_id],
# 		BD <- bd + bd_yr[yr_id],
# 		BS <- bs + bs_yr[yr_id],
# 		BG <- bg + bg_yr[yr_id],
# 		BA <- ba + ba_yr[yr_id],

# 		#adaptive priors
# 		a_cu[cu_id] ~ dnorm(0, sigma_cu),
# 		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

# 		#fixed priors
# 		a ~ dnorm(15,40),
# 		c(bd, bg, bs) ~ dnorm(0,1),
# 		ba ~ dnorm(15,40),
# 		sigma ~ dcauchy(0,2),
# 		sigma_cu ~ dcauchy(0,2),
# 		sigma_yr ~ dcauchy(0,2),
# 		Rho_yr ~ dlkjcorr(4)
# 		),
# 	data=d, 
# 	iter=5000, warmup=1000, chains=4, cores=4)

### Alternative formulationg with noncentered priors; more efficient but 
## difficult to interpret
# m1.dgs2b <- map2stan(
# 	alist(
# 		#likelihood
# 		TotalCount ~ dnorm(mu, sigma),
		
# 		#linear models
# 		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id,
# 		A <- a + a_cu[cu_id] + a_yr[yr_id],
# 		BD <- bd + bd_yr[yr_id],
# 		BS <- bs + bs_yr[yr_id],
# 		BG <- bg + bg_yr[yr_id],
# 		BA <- ba + ba_yr[yr_id],

# 		#adaptive priors
# 		a_cu[cu_id] ~ dnorm(0, sigma_cu),
# 		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr)[yr_id] ~ dmvnormNC(sigma_yr, Rho_yr),

# 		#fixed priors
# 		a ~ dnorm(15,50),
# 		c(bd, bg, bs) ~ dnorm(0,1),
# 		ba ~ dnorm(15,50),
# 		sigma ~ dcauchy(0,2),
# 		sigma_cu ~ dcauchy(0,2),
# 		sigma_yr ~ dcauchy(0,2),
# 		Rho_yr ~ dlkjcorr(2)
# 		),
# 	data=d, 
# 	iter=4000, warmup=1000, chains=4, cores=4)

# m1.dgs3 <- map2stan(
# 	alist(
# 		#likelihood
# 		TotalCount ~ dnorm(mu, sigma),
		
# 		#linear models
# 		mu <- a + a_cu[cu_id] + a_yr[yr_id] + bd*date_t + bg*growth_t + bs*fl_t 
# 		+ ba*age_id,
		
# 		#adaptive priors
# 		a_cu[cu_id] ~ dnorm(0, sigma_cu),
# 		a_yr[yr_id] ~ dnorm(0, sigma_yr),
		
# 		#fixed priors
# 		a ~ dnorm(15,50),
# 		c(bd, bg, bs) ~ dnorm(0,1),
# 		ba ~ dnorm(15, 50),
# 		sigma ~ dcauchy(0,2),
# 		sigma_cu ~ dcauchy(0,2),
# 		sigma_yr ~ dcauchy(0,2)
# 		),
# 	data=d, 
# 	iter=4000, warmup=1000, chains=4, cores=4)



################# Predictions for no random effects and one random slope models
# ### Link function for intercept and 1 slope model
# link.m1.dgs2 <- link(m1.dgs2, n=1000, data=d.pred,
# 	replace=list(a_cu=a_cu_sims, 
# 	# bg_cu=bg_cu_sims, bs_cu=bs_cu_sims, bd_cu=bd_cu_sims,
# 	a_yr=a_yr_sims, bg_yr=bg_yr_sims, bs_yr=bs_yr_sims, bd_yr=bd_yr_sims
# 	))

### Link function for intercept only model
# link.m1.dgs0 <- link(m1.dgs3, n=1000, data=d.pred,
# 	replace=list(a_cu=a_cu_sims, a_yr=a_yr_sims 
# 	))


## Posterior median and PI for intercept only
# mu.median <- apply(link.m1.dgs0, 2, median)
# lines(d.pred$date_t, mu.median, lwd=3)
# mu.PI <- apply(link.m1.dgs0, 2, PI, prob=0.9)
# shade(mu.PI, d.pred$date_t)

## Posterior median and PI for intercept and both slopes 
# mu.median <- apply(link.m1.dgs1$mu, 2, median)
# lines(d.pred$date_t, mu.median, lwd=1)
# mu.PI <- apply(link.m1.dgs1$mu, 2, PI, prob=0.9)
# shade(mu.PI, d.pred$date_t)





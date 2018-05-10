###### Migration Phenology - Simulation Sensitivity Analysis
###### Generate data with posterior estimates of fixed effects
###### and variance terms from fitted hierarchical models
###### Note: uses imported csv files of parameters listed on github 
###### estimates but could be readily modified to run with fresh runs
###### of models 
###### from supplementary data
###### February 2018
##################################################################
##################################################################


require(rethinking); require(MASS); require(dplyr); require(gplots); require(here)


dateModPars <- read.csv(here::here("JohnstoneStrait/dateModPars.csv"), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
countModPars <- read.csv(here::here("JohnstoneStrait/countModPars.csv"), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))

# ___________________________________________________________________
# Simulation sensitivity analysis 

## Generate data
## Define parameters based on estimates from models
## Import estimated fixed effects
a <- dateModPars[dateModPars$parameter=="alpha_fixed",2]
bs <- dateModPars[dateModPars$parameter=="betaSize_fixed",2]
br <- dateModPars[dateModPars$parameter=="betaRiver_fixed",2]
ba <- dateModPars[dateModPars$parameter=="betaAge0_fixed",2]
bb <- dateModPars[dateModPars$parameter=="betaAge2_fixed",2]
# CU covariance matrix
phiCUs <- dateModPars[c(58:61),2]
rhoCUs <- matrix(dateModPars[67:82,2], ncol=4, nrow=4)
covCU <- diag(phiCUs) %*% rhoCUs %*% diag(phiCUs)
# year covariance matrix
phiYrs <- dateModPars[c(62:66),2]
rhoYrs <- matrix(dateModPars[83:107,2], ncol=5, nrow=5)
covYr <- diag(phiYrs) %*% rhoYrs %*% diag(phiYrs)

## Random effects
### Simulate across 
nFishVec <- c(10, 50, 75)
nYrsVec <- c(3, 9)

dateSimPrecis <- vector("list", length(nFishVec))
datePrecis <- NULL

for(j in seq_along(nYrsVec)){ #loop through different year sample sizes
	set.seed(3333)
	# Define sample size
	nYrs <- nYrsVec[j]
	nCUs <- 9

	cuEffects <- mvrnorm(nCUs, rep(0, length.out=4), covCU)
	a_cu <- cuEffects[,1]
	bs_cu <- cuEffects[,2]
	ba_cu <- cuEffects[,3]
	bb_cu <- cuEffects[,4]
	
	yrEffects <- mvrnorm(nYrs, rep(0, length.out=5), covYr)
	a_yr <- yrEffects[,1]
	bs_yr <- yrEffects[,2]
	br_yr <- yrEffects[,3]
	ba_yr <- yrEffects[,4]
	bb_yr <- yrEffects[,5]
	
	
	# Simulate individuals within CUs and years
	#vectors of continuous variables
	yr_id <- rep(1:nYrs, each=(99999/9), length.out=99999)
	cu_id <- rep(1:nCUs, length.out=99999)
	fl_t <- rnorm(n=99999, mean=0, sd=1)
	river_t <- rnorm(n=99999, mean=0, sd=1)
	#vectors of categorical variables 
	age_id <- sample(c(0,1), 99999, replace=TRUE, prob=c(0.98,0.02)) #represents observed proportion of age 0
	age2_id <- age_id 
	for(i in seq_along(age2_id)){
		if(age2_id[i] == "0"){
			age2_id[i] <- 0
			draw <- runif(1, min=0, max=1)
			if(draw < 0.07){ #obs proportion of age 2 after removing age 0 
				age2_id[i] <- 1
			}
		} else{ age2_id[i] <- 0}
	}
	
	A <- a + a_cu[cu_id] + a_yr[yr_id]
	BS <- bs + bs_cu[cu_id] + bs_yr[yr_id]
	BR <- br + br_yr[yr_id]
	BA <- ba + ba_cu[cu_id] + ba_yr[yr_id]
	BB <- bb + bb_cu[cu_id] + bb_yr[yr_id]
	mu <- A +  BR*river_t + BS*fl_t + BA*age_id + BB*age2_id
			
	# Fixed sigma from model
	sigma <- dateModPars[6,2]
	
	# Generate data 
	entryDate <- rnorm(99999, mu, sigma)
	dateSimData <- data.frame(cu_id=cu_id, yr_id=yr_id, fl_t=fl_t, 
		river_t=river_t, age_id=age_id, age2_id=age2_id, entryDate=entryDate)
	dateSimData <- arrange(dateSimData, yr_id)

	for(i in 1:3){ # loop through different fish sample sizes
		set.seed(3333)
		nFish <- nFishVec[i]
		
		#subset simDat to appropriate length of years
		modelSimDat <- sample_n(dateSimData, (nFish*nYrs*nCUs))

		modRun <- map2stan(
			alist(
				#likelihood
				entryDate ~ dnorm(mu, sigma),
				
				#linear models
				mu <- A +  BR*river_t + BS*fl_t + BA*age_id + BB*age2_id,
				A <- a + a_cu[cu_id] + a_yr[yr_id],
				BR <- br + br_yr[yr_id],
				BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
				BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
				BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],
		
				#adaptive priors; internal model runs with non-centered priors to increase stability, but results qualitatively identical w/ standard
				c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnormNC(phi_cu, rho_cu), 
				c(a_yr, bs_yr, br_yr, ba_yr, bb_yr)[yr_id] ~ dmvnormNC(phi_yr, rho_yr),
		
				#fixed priors
				a ~ dnorm(120,30),
				ba ~ dnorm(5, 30),
				bb ~ dnorm(0, 30),
				c(bs,br) ~ dnorm(0,3),
				sigma ~ dcauchy(0,2),
				phi_cu ~ dcauchy(0,2),
				phi_yr ~ dcauchy(0,2),
				rho_cu ~ dlkjcorr(4),
				rho_yr ~ dlkjcorr(4)
			),
			data=modelSimDat, 
			iter=4500, warmup=750, chains=3, cores=3)  
	
		dateSimPrecis[[i]] <- precis(modRun, depth=2)
		names(dateSimPrecis)[[i]] <- paste(nFish, nYrs, sep="_")
		}
	datePrecis <- c(datePrecis, dateSimPrecis)
}


#### Plot results 
est3yr <- seq(from=93, to=107, by=1) #correspond to index slots in precis object equivalent to: "a","bs", "br", "ba", "bb", "sigma"
est9yr <- seq(from=123, to=137, by=1) #correspond to index slots in precis object equivalent to: "a","bs", "br", "ba", "bb","sigma"

parameter <- NULL
meanEffect <- NULL
lowerInterval <- NULL
upperInterval <- NULL
modRun <- NULL

for(i in unique(est3yr)){
	for(j in 1:3){
		parameter <- c(parameter, rownames(datePrecis[[j]]@output[i,]))
		meanEffect <- c(meanEffect,datePrecis[[j]]@output[i,1])
		lowerInterval <- c(lowerInterval, datePrecis[[j]]@output[i,3])
		upperInterval <- c(upperInterval, datePrecis[[j]]@output[i,4])
		modRun <- c(modRun, names(datePrecis)[[j]])
	}
}	
for(i in unique(est9yr)){	
	for(j in 4:6){
		parameter <- c(parameter, rownames(datePrecis[[j]]@output[i,]))
		meanEffect <- c(meanEffect,datePrecis[[j]]@output[i,1])
		lowerInterval <- c(lowerInterval, datePrecis[[j]]@output[i,3])
		upperInterval <- c(upperInterval, datePrecis[[j]]@output[i,4])
		modRun <- c(modRun, names(datePrecis)[[j]])
	}
}

# Clean up data frame prior to plotting
modelPars <- data.frame(par = parameter, mean = meanEffect, lowInt = lowerInterval, 
	upInt = upperInterval, modRun = modRun)
modelPars$par <- as.factor(modelPars$par)
modelPars$par = factor(modelPars$par,levels(modelPars$par)[c(1,5,4,2,3,6:15)])

# Index variables
trueEffects <- c(a,bs,br,ba,bb, sigma, phiCUs, phiYrs)
parSeq <- c("Intercept", "Beta Size", "Beta River", "Beta Age-1", "Beta Age-2", "Sigma",
	"Phi CU - Intercept", "Phi CU - Size", "Phi CU - Age-1", "Phi CU - Age-2",
	"Phi Yr - Intercept", "Phi Yr - Size", "Phi Yr - River", "Phi Yr - Age-1", "Phi Yr - Age-2"
	)
plotSeq <- levels(modelPars$par)

png("entryDateSim.png", width=6.5, height=6.5, units="in", res=400)
par(mar=c(2,3,2,0)+0.1, oma=c(1,0,0,0)+0.2, mfrow=c(4,4),  mgp=c(3, .5, 0))
for(i in seq_along(trueEffects)){
	dat <- modelPars[modelPars$par==plotSeq[i],]
	ylimits <- c(0.9*min(dat$lowInt),1.1*min(dat$lowInt))
	plotCI(c(0.9,1.9,2.9,1.1,2.1,3.1), dat$mean, ui=dat$upInt, li=dat$lowInt, xlim=c(0.5,3.5),
		pch=c(rep(21,length.out=3), rep(23,length.out=3)), main=parSeq[i], 
		ylab="", xlab="", xaxt='n',cex=1.2)
	abline(h=trueEffects[i])
	axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, 
		labels=c("","10","50","75",""), tcl=0) 
	if(i == 1 | i==5 | i==9 | i==13){
		mtext(side=2, line=1.6, text=expression(paste("Mean Estimate" )), 
			outer=FALSE, cex=0.8)
	}
	if(i == 12 | i==13 | i==14 | i==15){
		mtext(side=1, line=1.6, text=expression(paste("Mean Fish Per CU Yr" )), 
			outer=FALSE, cex=0.8)
	}
}
dev.off()


#____________________________________________________________________
### Part 2 - Simulate data for count model
a <- countModPars[countModPars$parameter=="alpha_fixed",2]
bd <- countModPars[countModPars$parameter=="betaDate_fixed",2]
bs <- countModPars[countModPars$parameter=="betaSize_fixed",2]
bg <- countModPars[countModPars$parameter=="betaGrowth_fixed",2]
ba <- countModPars[countModPars$parameter=="betaAge0_fixed",2]
bb <- countModPars[countModPars$parameter=="betaAge2_fixed",2]
# CU covariance matrix
phiCUs <- countModPars[c(80:85),2]
rhoCUs <- matrix(countModPars[92:127,2], ncol=6, nrow=6)
covCU <- diag(phiCUs) %*% rhoCUs %*% diag(phiCUs)
# year covariance matrix
phiYrs <- countModPars[c(86:91),2]
rhoYrs <- matrix(countModPars[128:163,2], ncol=6, nrow=6)
covYr <- diag(phiYrs) %*% rhoYrs %*% diag(phiYrs)



### Simulate with different sample sizes 
nFishVec <- c(10, 50, 75)
nYrsVec <- c(3, 9)

simPrecis <- vector("list", length(nFishVec))
countPrecis <- NULL

for(j in seq_along(nYrsVec)){
	set.seed(3333)
	# Define sample size
	nYrs <- nYrsVec[j]
	nCUs <- 9
	
	# Generate dataset based on how many years of effort
	## Random effects
	cuEffects <- mvrnorm(nCUs, rep(0, length.out=6), covCU)
	a_cu <- cuEffects[,1]
	bd_cu <- cuEffects[,2]
	bs_cu <- cuEffects[,3]
	bg_cu <- cuEffects[,4]
	ba_cu <- cuEffects[,5]
	bb_cu <- cuEffects[,6]
	
	yrEffects <- mvrnorm(nYrs, rep(0, length.out=6), covYr)
	a_yr <- yrEffects[,1]
	bd_yr <- yrEffects[,2]
	bs_yr <- yrEffects[,3]
	bg_yr <- yrEffects[,4]
	ba_yr <- yrEffects[,5]
	bb_yr <- yrEffects[,6]
	
	# Simulate individuals within CUs and years
	#vectors of continuous variables
	yr_id <- rep(1:nYrs, each=(99999/9), length.out=99999)
	cu_id <- rep(1:nCUs, length.out=99999)
	date_t <- rnorm(n=99999, mean=0, sd=1)
	fl_t <- rnorm(n=99999, mean=0, sd=1)
	growth_t <- rnorm(n=99999, mean=0, sd=1)
	#vectors of categorical variables
	age_id <- sample(c(0,1), 99999, replace=TRUE, prob=c(0.98,0.02)) 
	age2_id <- age_id 
	for(i in seq_along(age2_id)){
		if(age2_id[i] == "0"){
			age2_id[i] <- 0
			draw <- runif(1, min=0, max=1)
			if(draw < 0.07){
				age2_id[i] <- 1
			}
		} else{ age2_id[i] <- 0}
	}
	
	A <- a + a_cu[cu_id] + a_yr[yr_id]
	BD <- bd + bd_cu[cu_id] + bd_yr[yr_id]
	BS <- bs + bs_cu[cu_id] + bs_yr[yr_id]
	BG <- bg + bg_cu[cu_id] + bg_yr[yr_id]
	BA <- ba + ba_cu[cu_id] + ba_yr[yr_id]
	BB <- bb + bb_cu[cu_id] + bb_yr[yr_id]
	mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id + BB*age2_id
	
	# Fixed sigma from model
	sigma <- countModPars[7,2]
	
	# Generate data 
	migDuration <- rnorm(99999, mu, sigma)
	simData <- data.frame(cu_id=cu_id, yr_id=yr_id, date_t=date_t, fl_t=fl_t, 
		growth_t=growth_t, age_id=age_id, age2_id=age2_id, migDuration=migDuration)
	
	# Sample based on 
	for(i in seq_along(nFishVec)){
		nFish <- nFishVec[i]
		
		#subset simDat to appropriate length of years
		set.seed(3333)
		modelSimDat <- sample_n(simData, (nFish*nYrs*nCUs))

		modRun <- map2stan(
			alist(
				#likelihood
				migDuration ~ dnorm(mu, sigma),
				#linear models
				mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id + BB*age2_id,
				A <- a + a_cu[cu_id] + a_yr[yr_id],
				BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
				BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
				BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
				BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
				BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],
				#adaptive priors; internal model runs with non-centered priors to increase stability, but results qualitatively identical w/ standard
				c(a_cu, bd_cu, bs_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnormNC(phi_cu, rho_cu),
				c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnormNC(phi_yr, rho_yr),
				#fixed priors
				a ~ dnorm(15,30),
				c(bd, bg, bs) ~ dnorm(0,3),
				ba ~ dnorm(5, 30),
				bb ~ dnorm(0, 30),
				sigma ~ dcauchy(0,2),
				phi_cu ~ dcauchy(0,2),
				phi_yr ~ dcauchy(0,2),
				rho_cu ~ dlkjcorr(4),
				rho_yr ~ dlkjcorr(4)
				),
			data=modelSimDat, 
			iter=4000, warmup=750, chains=3, cores=3)  
	
		simPrecis[[i]] <- precis(modRun, depth=2)
		names(simPrecis)[[i]] <- paste(nFish, nYrs, sep="_")
	}
	countPrecis <- c(countPrecis, simPrecis)
}

est3yr <- seq(from=145, to=163, by=1) #correspond to index slots in precis object equivalent to: "a","bd","bs","bg","ba","bb","sigma"
est9yr <- seq(from=181, to=199, by=1)

parameter <- NULL
meanEffect <- NULL
lowerInterval <- NULL
upperInterval <- NULL
modRun <- NULL

for(i in unique(est3yr)){
	for(j in 1:3){
		parameter <- c(parameter, rownames(countPrecis[[j]]@output[i,]))
		meanEffect <- c(meanEffect,countPrecis[[j]]@output[i,1])
		lowerInterval <- c(lowerInterval, countPrecis[[j]]@output[i,3])
		upperInterval <- c(upperInterval, countPrecis[[j]]@output[i,4])
		modRun <- c(modRun, names(countPrecis)[[j]])
	}
}	
for(i in unique(est9yr)){	
	for(j in 4:6){
		parameter <- c(parameter, rownames(countPrecis[[j]]@output[i,]))
		meanEffect <- c(meanEffect,countPrecis[[j]]@output[i,1])
		lowerInterval <- c(lowerInterval, countPrecis[[j]]@output[i,3])
		upperInterval <- c(upperInterval, countPrecis[[j]]@output[i,4])
		modRun <- c(modRun, names(countPrecis)[[j]])
	}
}

# Clean up data frame prior to plotting
modelPars <- data.frame(par = parameter, mean = meanEffect, lowInt = lowerInterval, 
	upInt = upperInterval, modRun = modRun)
modelPars$par <- as.factor(modelPars$par)
modelPars$yr <- 1
for(i in 1:nrow(modelPars)){
	if(modelPars$modRun[i] == "10_3" | modelPars$modRun[i] == "50_3" | modelPars$modRun[i] == "75_3"){
		modelPars$yr[i] <- 0
	}
}
modelPars$par = factor(modelPars$par,levels(modelPars$par)[c(1,4,6,5,2,3,7:19)])

# Index variables
trueEffects <- c(c(a,bd,bs,bg,ba,bb,sigma), phiCUs, phiYrs)
parSeq <- c("Intercept", "Beta Date", "Beta Size", "Beta Growth", "Beta Age-1", "Beta Age-2", "Sigma",
	"Phi CU - Intercept", "Phi CU - Date", "Phi CU - Size", "Phi CU - Growth", "Phi CU - Age-1", "Phi CU - Age-2",
	"Phi Yr - Intercept", "Phi Yr - Date", "Phi Yr - Size", "Phi Yr - Growth", "Phi Yr - Age-1", "Phi Yr - Age-2"
	)
plotSeq <- levels(modelPars$par)

png("countSim.png", width=6.5, height=6.5, units="in", res=400)
par(mar=c(2,3,2,0)+0.1, oma=c(1,0,0,0)+0.2, mfrow=c(5,4),  mgp=c(3, .5, 0))
for(i in seq_along(trueEffects)){
	dat <- modelPars[modelPars$par==plotSeq[i],]
	ylimits <- c(0.9*min(dat$lowInt),1.1*min(dat$lowInt))
	plotCI(c(0.9,1.9,2.9,1.1,2.1,3.1), dat$mean, ui=dat$upInt, li=dat$lowInt, xlim=c(0.5,3.5),
		pch=c(rep(21,length.out=3), rep(23,length.out=3)), main=parSeq[i], 
		ylab="", xlab="", xaxt='n',cex=1.2)
	abline(h=trueEffects[i])
	axis(1, tick=T, at=c(0,1,2,3,4), lwd=2, lwd.ticks=1, cex.axis=1, 
		labels=c("","10","50","75",""), tcl=0) 
	if(i == 1 | i==5 | i==9 | i==13 | i==17){
		mtext(side=2, line=1.6, text=expression(paste("Mean Estimate" )), 
			outer=FALSE, cex=0.8)
	}
	if(i == 16 | i==17 | i==18 | i==19){
		mtext(side=1, line=1.6, text=expression(paste("Mean Fish Per CU Yr" )), 
			outer=FALSE, cex=0.8)
	}
}
dev.off()

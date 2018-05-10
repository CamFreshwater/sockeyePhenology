##################################################################
##################################################################
######## Migration Phenology - Hierarchical Model 
###### Use top ranked models to explore individual, population and 
###### year scale effects on phenology timing and duration
###### February 2018
##################################################################
##################################################################

require(rethinking)
require(here) #package supporting generic import assuming your selected wd contains file of interest

phenDat <- read.csv(here("individualOtolithData.csv"), stringsAsFactors = FALSE, strip.white = TRUE, 
	na.strings = c("NA",""))


#____________________________________________________________________
## Data clean

phenDat$Age <- as.factor(phenDat$age)

#### Make index variables
phenDat$cu_id <- coerce_index(phenDat$cu)
phenDat$yr_id <- coerce_index(phenDat$year)
phenDat$age_id <- 0
phenDat$age_id[phenDat$Age=="0"] <- 1
phenDat$age2_id <- 0
phenDat$age2_id[phenDat$Age=="2"] <- 1

## Transform explanatory variables
phenDat$fl_t <- (phenDat$entryFL - mean(phenDat$entryFL))/sd(phenDat$entryFL) 
phenDat$date_t <- (phenDat$oceanEntryDate - mean(phenDat$oceanEntryDate))/sd(phenDat$oceanEntryDate)
phenDat$growth_t <- (phenDat$meanGrowth - mean(phenDat$meanGrowth))/sd(phenDat$meanGrowth)
phenDat$river_t <- (phenDat$riverDist - mean(phenDat$riverDist))/sd(phenDat$riverDist)

## Remove extra columns
d <- phenDat[,c("yr_id","cu_id","age_id","age2_id","TotalCount","EntryDate","fl_t","date_t",
	"growth_t","river_t")]



#____________________________________________________________________
## Models

## Total count model
modCount <- map2stan(
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
		c(a_cu, bd_cu, bs_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),
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
	data=d, 
	iter=9000, warmup=1000, chains=3, cores=3)


### Entry date models
modDate <- map2stan(
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
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bs_yr, br_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr,	rho_yr),

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
	data=d, 
	iter=9000, warmup=1000, chains=3, cores=3)


### Check convergence
precis(modCount, depth=2)
plot(modCount)
precis(modDate, depth=2)
plot(modDate)

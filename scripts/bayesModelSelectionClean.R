##################################################################
##################################################################
######## Migration Phenology - Hierarchical Model 
###### Build and select top Bayesian hierarchical models 
###### to estimate effects of entry date, entry size, and growth rate 
###### on sockeye phenology (days of residence and midday) while 
###### accounting for interannual and population specific effects;
###### EDIT 1 Aug 14, 2017 to include age 2 effects
###### EDIT 2 Oct 19, 2017 to weaken priors
##################################################################
##################################################################

require(rethinking)
require(car)

phenDat <- read.csv(here::here("individualOtolithData.csv"), stringsAsFactors = FALSE, strip.white = TRUE, 
	na.strings = c("NA",""))


####### DATA CLEAN #######

phenDat$Age <- as.factor(phenDat$Age)

## Make index variables
phenDat$cu_id <- coerce_index(phenDat$CU)
phenDat$yr_id <- coerce_index(phenDat$Year)
phenDat$age_id <- 0
phenDat$age_id[phenDat$Age=="0"] <- 1
phenDat$age2_id <- 0
phenDat$age2_id[phenDat$Age=="2"] <- 1

## Transform explanatory variables
phenDat$fl_t <- (phenDat$EntryFL - mean(phenDat$EntryFL))/sd(phenDat$EntryFL) 
phenDat$date_t <- (phenDat$EntryDate - mean(phenDat$EntryDate))/sd(phenDat$EntryDate)
phenDat$growth_t <- (phenDat$BodyGrowth - mean(phenDat$BodyGrowth))/sd(phenDat$BodyGrowth)
phenDat$river_t <- (phenDat$RiverDist - mean(phenDat$RiverDist))/sd(phenDat$RiverDist)
phenDat$lakelat_t <- (phenDat$LakeLat - mean(phenDat$LakeLat))/sd(phenDat$LakeLat)


## Remove extra columns
d <- phenDat[,c("FishNumber","yr_id","cu_id","age_id","age2_id","TotalCount","EntryDate",
	"fl_t","date_t","growth_t","river_t","lakelat_t")]


vif(lm(TotalCount ~ age2_id + age_id + fl_t + date_t + growth_t, data=d))
vif(lm(EntryDate ~ age2_id + age_id + fl_t + river_t, data=d))

################################################################


############## BUILD MODELS ##############


######## Migration duration models
## No covariate model
m1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)
		
## Single covariate models
m1.d1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bd_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		bd ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m1.g1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BG*growth_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		bg ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m1.s1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		bs ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

### Two covariate model 
m1.dg1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bd_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		c(bd, bg) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m1.ds1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bd_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		c(bd, bs) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m1.gs1 <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BG*growth_t + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bg_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bg_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		c(bg, bs) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

### Three covariate model - ranked as best
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
		c(a_cu, bd_cu, bs_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		c(bd, bg, bs) ~ dnorm(0,3),
		c(ba,bb) ~ dnorm(15, 30),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

###### Check model convergence
## n_eff and Rhat
precis(m1, depth=2)
precis(m1.d1, depth=2)
precis(m1.g1, depth=2) ## some issues with a estimates until priors inc
precis(m1.s1, depth=2)
precis(m1.dg1, depth=2)
precis(m1.ds1, depth=2)
precis(m1.gs1, depth=2)
precis(m1.dgs1, depth=2)

###### Model comparison
compare(m1.d1, m1.s1, m1.g1, m1.dg1, m1.ds1, m1.gs1, m1.dgs1, m1)
compare(m1, m1.dgs1)

## Saturated models have best fit (difference = 13, 100% of weight)

plot(coeftab(m1.d, m1.s, m1.g, m1.dg, m1.ds, m1.gs, m1.dgs))
## Parameter estimates similar across models

######################################################
######################################################



######## Entry date models

## Single covariate models
m3 <- map2stan(
	alist(
		#likelihood
		EntryDate ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m3.s1 <- map2stan(
	alist(
		#likelihood
		EntryDate ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		bs ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m3.r1 <- map2stan(
	alist(
		#likelihood
		EntryDate ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BR*river_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BR <- br + br_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, br_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		br ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

m3.l1 <- map2stan(
	alist(
		#likelihood
		EntryDate ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BL*lakelat_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BL <- bl + bl_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bl_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		bl ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

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
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bs_yr, br_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

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
	iter=8000, warmup=1000, chains=3, cores=3)

m3.ls1 <- map2stan(
	alist(
		#likelihood
		EntryDate ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A +  BL*lakelat_t + BS*fl_t + BA*age_id + BB*age2_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BL <- bl + bl_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BA <- ba + ba_cu[cu_id] + ba_yr[yr_id],
		BB <- bb + bb_cu[cu_id] + bb_yr[yr_id],

		#adaptive priors
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, phi_cu, rho_cu),
		c(a_yr, bs_yr, bl_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, phi_yr, rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(5, 30),
		bb ~ dnorm(0, 30),
		c(bs,bl) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		phi_cu ~ dcauchy(0,2),
		phi_yr ~ dcauchy(0,2),
		rho_cu ~ dlkjcorr(4),
		rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=8000, warmup=1000, chains=3, cores=3)

###### Check model convergence
## n_eff and Rhat
precis(m3, depth=2)
precis(m3.s1, depth=2)
precis(m3.r1, depth=2)
precis(m3.rs1, depth=2)

###### Model comparison
compare(m3, m3.s1, m3.r1, m3.l1, m3.rs1, m3.ls1)



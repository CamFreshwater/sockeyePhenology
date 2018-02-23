##################################################################
##################################################################
######## Johnstone Strait - Hierarchical Modeling Comparison
###### Build and select top Bayesian hierarchical models 
###### to estimate effects of entry date, entry size, and growth rate 
###### on sockeye phenology (days of residence and midday) while 
###### accounting for interannual and population specific effects;
###### edited again Aug 14 to include age 2 effects
###### June 2017
##################################################################
##################################################################


rm(list=ls())


setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")



require(rethinking)
require(car)


JSdat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))



####### DATA CLEAN #######

JSdat$Age <- as.factor(JSdat$Age)

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
JSdat$lakelat_t <- (JSdat$LakeLat - mean(JSdat$LakeLat))/sd(JSdat$LakeLat)


## Remove extra columns
d <- JSdat[,c("FishNumber","yr_id","cu_id","age_id","age2_id","TotalCount","EntryDate",
	"fl_t","date_t","growth_t","river_t","lakelat_t")]


# plot(TotalCount ~ fl_t, data=d)
vif(lm(TotalCount ~ age2_id + age_id + fl_t + date_t + growth_t, data=d))
vif(lm(EntryDate ~ age2_id + age_id + fl_t + river_t, data=d))

################################################################
################################################################
################################################################


############## BUILD MODELS ##############



######## Total count models

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
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=6000, warmup=1000, chains=3, cores=3)
		
	
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
		c(a_cu, bd_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		bd ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=6000, warmup=1000, chains=3, cores=3)


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
		c(a_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		bg ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=6000, warmup=1000, chains=3, cores=3)



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
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		bs ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)



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
		c(a_cu, bd_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		c(bd, bg) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)



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
		c(a_cu, bd_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		c(bd, bs) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)



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
		c(a_cu, bg_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bg_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		c(bg, bs) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)




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
		c(a_cu, bd_cu, bs_cu, bg_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		c(bd, bg, bs) ~ dnorm(0,1),
		c(ba,bb) ~ dnorm(15, 30),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7500, warmup=1000, chains=3, cores=3)


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
# precis(m1.d2, depth=2)
# precis(m1.g2, depth=2) ## some issues with a estimates until priors inc
# precis(m1.s2, depth=2)
# precis(m1.dg2, depth=2)
# precis(m1.ds2, depth=2)
# precis(m1.gs2, depth=2)
# precis(m1.dgs2, depth=2)
# precis(m1.dgs3, depth=2)



###### Model comparison
compare(m1.d1, m1.s1, m1.g1, m1.dg1, m1.ds1, m1.gs1, m1.dgs1)
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
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=5000, warmup=1000, chains=3, cores=3)


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
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bs_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		bs ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=5000, warmup=1000, chains=3, cores=3)


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
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, br_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		br ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)

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
		c(a_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bl_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		bl ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)


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
		c(bs,br) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)

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
		c(a_cu, bs_cu, ba_cu, bb_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bs_yr, bl_yr, ba_yr, bb_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(120,30),
		ba ~ dnorm(15, 30),
		bb ~ dnorm(0, 30),
		c(bs,bl) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=7000, warmup=1000, chains=3, cores=3)



###### Check model convergence
## n_eff and Rhat
precis(m3, depth=2)
precis(m3.s1, depth=2)
precis(m3.r1, depth=2)
precis(m3.rs1, depth=2)
# precis(m3.rs2, depth=2)

###### Model comparison
compare(m3, m3.s1, m3.r1, m3.l1, m3.rs1, m3.ls1)



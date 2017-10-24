##################################################################
##################################################################
######## Johnstone Strait - Preliminary Modeling 
###### Compare frequentist and Bayesian multilevel, hierarchical 
###### models to examine changes and ID best distribution
###### June 2017
##################################################################
##################################################################



setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait")



require(rethinking)
require(car)
require(nlme)
require(lme4)


JSdat <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatCLEAN.csv", 
          stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))



####### DATA CLEAN #######

JSdat$Age <- as.factor(JSdat$Age)

## Make index variables
JSdat$cu_id <- coerce_index(JSdat$CU)
JSdat$yr_id <- coerce_index(JSdat$Year)
JSdat$age_id <- 0
JSdat$age_id[JSdat$Age=="0"] <- 1

## Transform explanatory variables
JSdat$fl_t <- (JSdat$EntryFL - mean(JSdat$EntryFL))/sd(JSdat$EntryFL) 
JSdat$date_t <- (JSdat$EntryDate - mean(JSdat$EntryDate))/sd(JSdat$EntryDate)
JSdat$growth_t <- (JSdat$BodyGrowth - mean(JSdat$BodyGrowth))/sd(JSdat$BodyGrowth)


## Remove extra columns
d <- JSdat[,c("FishNumber","yr_id","cu_id","age_id","TotalCount","EntryDate",
	"fl_t","date_t","growth_t")]


## Round total count variable
d$TotalCount <- round(d$TotalCount, digits=0) 
d$EntryDate <- round(d$EntryDate, digits=0) 


plot(TotalCount ~ growth_t, data=d, col=as.numeric(Age), pch=16)


##################################################


#### Check for collinearity using basic frequentist model
vif(lm(TotalCount ~ fl_t + date_t + growth_t + age_id, data=d))
## Good to go


####### BUILD MODELS #######

#### Practice frequentist models and compare to Bayesian equivalents
#### using single covariate

plot( ~ date_t, data=d)


##### Example 1. Assumed Gaussian distribution

### Frequentist
mod <- lm(TotalCount ~ date_t, data=d)

summary(mod)
hist(resid(mod))
plot(mod)

## Bayesian equivalent w/ Gaussian
m1.prac <- map2stan(
	alist(
		TotalCount ~ dnorm(mu, sigma),
		mu ~ a + bD*date_t,
		a ~ dnorm(0,100),
		bD ~ dnorm(0,1),
		sigma ~ dcauchy(0,2)
		),
	data=d, 
	iter=3000, warmup=1000, chains=3, cores=3)

plot(precis(m1.prac))
precis(m1.prac)
##############################################



##### Example 2. Assumed Poisson distribution

### Frequentist
mod2 <- glm(TotalCount ~ date_t, family=poisson, data=d)

summary(mod2)
hist(resid(mod2))
plot(mod2)


## Bayesian equivalent w/ Gaussian
m2.prac <- map2stan(
	alist(
		TotalCount ~ dpois(lambda),
		log(lambda) ~ a + bD*date_t,
		a ~ dnorm(0,10),
		bD ~ dnorm(0,1)
		),
	data=d, 
	iter=3000, warmup=1000, chains=3, cores=3)

plot(precis(m2.prac))
##############################################



##### Example 3. Assumed Poisson distribution w/ 
#### random intercepts

### Frequentist
mod3 <- glmer(TotalCount ~ date_t + (1|cu_id), family=poisson, data=d)
mod3b <- lmer(TotalCount ~ date_t + (1|cu_id), data=d)


summary(mod3)
summary(mod3b)
hist(resid(mod3))
hist(resid(mod3b))
plot(mod3b)

shapiro.test(resid(mod3))
shapiro.test(resid(mod3b))


## Bayesian equivalents w/ hierarchical effects
## Gamma poisson distribution
m3.prac <- map2stan(
	alist(
		TotalCount ~ dgampois(mu, scale),
		log(mu) ~ a + a_cu[cu_id] + bS*fl_t + bA*age_id,
		a_cu[cu_id] ~ dnorm(0, tau),
		a ~ dnorm(5,10),
		bA ~ dnorm(1,10),
		bS ~ dnorm(0,3),
		scale ~ dcauchy(0,2),
		tau ~ dcauchy(0,2)
		),
	data=d, 
	constraints=list(scale="lower=0"),
	iter=3000, warmup=1000, chains=3, cores=3)

plot(m3.prac)
precis(m3.prac, depth=2)
plot(precis(m3.prac, depth=2))


## Poisson distribution
m3b.prac <- map2stan(
	alist(
		TotalCount ~ dpois(lambda),
		log(lambda) ~ a + a_cu[cu_id] + bS*fl_t + bA*age_id,
		a_cu[cu_id] ~ dnorm(0, tau),
		a ~ dnorm(5,10),
		bA ~ dnorm(1,10),
		bS ~ dnorm(1,1),
		tau ~ dcauchy(0,2)
		),
	data=d, 
	iter=3000, warmup=1000, chains=3, cores=3)

plot(m3b.prac)
plot(precis(m3b.prac,depth=2))
precis(m3b.prac, depth=2)


## Gaussian distribution
m3c.prac <- map2stan(
	alist(
		TotalCount ~ dnorm(mu, sigma),
		mu ~ a + a_cu[cu_id] + bS*fl_t + bA*age_id,
		a_cu[cu_id] ~ dnorm(0, tau),
		a ~ dnorm(10,100),
		bA ~ dnorm(10,100),
		bS ~ dnorm(1,1),
		sigma ~ dcauchy(0,2),
		tau ~ dcauchy(0,2)
		),
	data=d, 
	iter=3000, warmup=1000, chains=3, cores=3)

plot(m3c.prac)
plot(precis(m3c.prac,depth=2))
precis(m3c.prac, depth=2)



##### Compare predicted values from simple Gaussian
#### gamma-poisson and poisson models

compare(m3.prac, m3b.prac, m3c.prac)
## Gaussian model has much superior fit


### Extract samples
postA <- extract.samples(m3.prac)
postB <- extract.samples(m3b.prac)
postC <- extract.samples(m3c.prac)


d.pred <- list(
	# explanatory variable of interest
	fl_t = seq(from=-4.5, to=4.5, length.out=50),
	age_id = rep(0, 50),
	# random effects
	cu_id = rep(1, 50)
	)

# ## Draw slopes from distribution using estimates of variation
# ## in cu and yr

## Gamma poisson
a_cu_sims <- matrix(rnorm(9000, 0, postA$tau))
link.m3a.prac <- link(m3.prac, n=1000, data=d.pred,
	replace=list(a_cu=a_cu_sims))

## Poisson
a_cu_sims <- matrix(rnorm(9000, 0, postA$tau))
link.m3b.prac <- link(m3b.prac, n=1000, data=d.pred,
	replace=list(a_cu=a_cu_sims))

## Gaussian
a_cu_sims <- matrix(rnorm(9000, 0, postA$tau))
link.m3c.prac <- link(m3c.prac, n=1000, data=d.pred,
	replace=list(a_cu=a_cu_sims))


par(mfrow=c(2,2))
## Raw data
plot(TotalCount ~ fl_t, data=d, pch=16, xlab="Entry Size",
	ylab="Days Resident")
## Posterior median for each distribuiton
mu.median <- apply(link.m3a.prac, 2, median)
lines(d.pred$fl_t, mu.median)
## Plot 90% CI
mu.PI <- apply(link.m3a.prac, 2, PI, prob=0.9)
shade(mu.PI, d.pred$fl_t)

## Raw data
plot(TotalCount ~ fl_t, data=d, pch=16, xlab="Entry Size",
	ylab="Days Resident")
mu.medianB <- apply(link.m3b.prac, 2, median)
lines(d.pred$fl_t, mu.medianB)
## Plot 90% CI
mu.PI <- apply(link.m3b.prac, 2, PI, prob=0.9)
shade(mu.PI, d.pred$fl_t)

## Raw data
plot(TotalCount ~ fl_t, data=d, pch=16, xlab="Entry Size",
	ylab="Days Resident")
mu.medianC <- apply(link.m3c.prac, 2, median)
lines(d.pred$fl_t, mu.medianC)
## Plot 90% CI
mu.PI <- apply(link.m3c.prac, 2, PI, prob=0.9)
shade(mu.PI, d.pred$fl_t)


### Gaussian has tighter confidence intervals and much better score
## no reason to use others




#### STAGE 1 CONCLUSIONS - frequentist and Bayesian models
#### provide similar estimates; m3e seems ideal for interpretability and 
#### consistency will use gaussian distribution for both

##############################################





##### Example 4. Assumed Gaussian distribution w/ 
#### random intercepts and random slopes


### Frequentist
mod4 <- lmer(TotalCount ~ date_t + (1 + date_t|CU), data=d)


summary(mod4)
hist(resid(mod4))
plot(mod4)

## Check residuals for practice
E2 <- resid(mod4, type="response")
F2 <- fitted(mod4)
plot(E2 ~ F2)
plot(E2 ~ d$date_t)

## Gaussian distribution
m4.prac <- map2stan(
	alist(
		TotalCount ~ dnorm(mu, sigma),
		mu ~ a_cu[cu_id] + bD_cu[cu_id]*date_t,
		c(a_cu,bD_cu)[cu_id] ~ dmvnorm2(c(a, bD), sigma_cu, Rho),
		a ~ dnorm(0,100),
		bD ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		Rho ~ dlkjcorr(2)
		),
	data=d, 
	iter=3000, warmup=1000, chains=4, cores=3)

plot(m4.prac)
plot(precis(m4.prac,depth=2))
precis(m4.prac, depth=2)

##############################################



##### Example 5. Assumed Gaussian distribution w/ 
#### NESTED, random intercepts 


### Frequentist
mod5 <- lmer(TotalCount ~ date_t + (1 |Year/CU), data=d)


summary(mod5)
hist(resid(mod5))
plot(mod5)


## Gaussian distribution 
## (NOTE THAT SINGLE GRAND MEAN IS ESTIMATED FOR INTERPRETABILITY)
m5.prac <- map2stan(
	alist(
		TotalCount ~ dnorm(mu, sigma),
		mu ~ a + a_cu[cu_id] + a_yr[yr_id] + bD*date_t,
		a_cu[cu_id] ~ dnorm(0, sigma_cu),
		a_yr[yr_id] ~ dnorm(0, sigma_yr),
		a ~ dnorm(0,100),
		bD ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2)
		),
	data=d, 
	iter=3000, warmup=1000, chains=4, cores=3)

plot(m5.prac)
plot(precis(m5.prac,depth=2))
precis(m5.prac, depth=2)

##############################################



##### Example 6. FINAL model structure w/ nested
#### random intercepts and slopes for both yr and CU  


### Frequentist
mod6 <- lmer(TotalCount ~ fl_t + date_t + growth_t +
	 (1 + date_t + fl_t + growth_t|yr_id/cu_id), data=d, REML=TRUE, family=poisson)
summary(mod6)
mod6b <- lmer(TotalCount ~ date_t + growth_t +
	 (1 + date_t + growth_t |yr_id/cu_id), data=d, REML=TRUE)

AIC(mod6, mod6b)

summary(mod6)
hist(resid(mod6))
plot(mod6)


mod5 <- lmer(EntryDate ~ fl_t + age_id + (1 + age_id + fl_t|cu_id), data=d)

summary(mod5)
summary(mod6)
		


## Gaussian distribution 
## NOTE: SINGLE GRAND MEAN IS ESTIMATED FOR INTERPRETABILITY
## NOTE2: Unique model structure for simplicity, separate out 
## linear models and priors
m6.prac <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bg_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bg_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		c(bd,bg) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=5000, warmup=1000, chains=3, cores=3)

plot(m6.prac)
plot(precis(m6.prac,depth=2))
precis(m6.prac, depth=2)
precis(m6.prac)

post <- extract.samples(m6.prac)
str(post)
pairs(post)



#### Try poisson distribution
## Poisson distribution
m6b.prac <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dpois(lambda),
		
		#linear models
		log(lambda) <- A + BD*date_t + BG*growth_t,
		
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bg_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bg_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,30),
		c(bd,bg) ~ dnorm(0,3),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,2),
		sigma_yr ~ dcauchy(0,2),
		Rho_cu ~ dlkjcorr(4),
		Rho_yr ~ dlkjcorr(4)
		),
	data=d, 
	iter=5000, warmup=1000, chains=3, cores=3)


#### Compare both models
compare(m6.prac, m6b.prac)

#### Plot predictions from gaussian and poisson model

		
##############################################

### Three covariate model - ranked as best
## Compare three versions: a) without age-0 parameter, b) with age-0
## c) with slopes only on years and d) with only intercepts


###### CONCLUSIONS - more support for one slope than both, but not full
##### include both model types in selection

#### Total count
 m1.dgsA <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bs_cu, bg_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,50),
		c(bs, bg, bd) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3),
		Rho_cu ~ dlkjcorr(3),
		Rho_yr ~ dlkjcorr(3)
		),
	data=d, 
	iter=4000, warmup=1000, chains=4, cores=4)


m1.dgsB <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
		BA <- ba + ba_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bs_cu, bg_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,50),
		c(bs, bg, bd) ~ dnorm(0,1),
		ba ~ dnorm(15, 50),
		sigma ~ dcauchy(0,3),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3),
		Rho_cu ~ dlkjcorr(3),
		Rho_yr ~ dlkjcorr(3)
		),
	data=d, 
	iter=4000, warmup=1000, chains=4, cores=4)

precis(m1.dgsA, depth=2)
precis(m1.dgsB, depth=2)


m1.dgs.noslope <- map2stan(
	alist(
		TotalCount ~ dnorm(mu, sigma),
		mu ~ a + a_cu[cu_id] + a_yr[yr_id] + bs*fl_t + bd*date_t 
			+ bg*growth_t + ba*age_id,
		a_cu[cu_id] ~ dnorm(0, sigma_cu),
		a_yr[yr_id] ~ dnorm(0, sigma_yr),
		a ~ dnorm(15,50),
		ba ~ dnorm(15, 50),
		bs ~ dnorm(0,1),
		bg ~ dnorm(0,1),
		bd ~ dnorm(0,1),
		sigma ~ dcauchy(0,3),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3)
		),
	data=d, 
	iter=3000, warmup=1000, chains=4, cores=3)


m1.dgs.oneslope <- map2stan(
	alist(
		#likelihood
		TotalCount ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_yr[yr_id],
		BS <- bs + bs_yr[yr_id],
		BG <- bg + bg_yr[yr_id],
		BA <- ba + ba_yr[yr_id],

		#adaptive priors
		a_cu[cu_id] ~ dnorm(0, sigma_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,50),
		c(bs, bg, bd) ~ dnorm(0,1),
		ba ~ dnorm(15, 50),
		sigma ~ dcauchy(0,3),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3),
		Rho_yr ~ dlkjcorr(3)
		),
	data=d, 
	iter=4000, warmup=1000, chains=4, cores=4)


compare(m1.dgsA, m1.dgsB, m1.dgs.noslope, m1.dgs.oneslope)

precis(m1.dgs.noslope,depth=2)
precis(m1.dgs.oneslope,depth=2)




### Midday
 m2.dgsA <- map2stan(
	alist(
		#likelihood
		MidDay ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bs_cu, bg_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,50),
		c(bs, bg, bd) ~ dnorm(0,1),
		sigma ~ dcauchy(0,2),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3),
		Rho_cu ~ dlkjcorr(3),
		Rho_yr ~ dlkjcorr(3)
		),
	data=d, 
	iter=4000, warmup=1000, chains=4, cores=4)


m2.dgsB <- map2stan(
	alist(
		#likelihood
		MidDay ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_cu[cu_id] + bd_yr[yr_id],
		BS <- bs + bs_cu[cu_id] + bs_yr[yr_id],
		BG <- bg + bg_cu[cu_id] + bg_yr[yr_id],
		BA <- ba + ba_yr[yr_id],

		#adaptive priors
		c(a_cu, bd_cu, bs_cu, bg_cu)[cu_id] ~ dmvnorm2(0, sigma_cu, Rho_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,50),
		c(bs, bg, bd) ~ dnorm(0,1),
		ba ~ dnorm(15, 50),
		sigma ~ dcauchy(0,3),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3),
		Rho_cu ~ dlkjcorr(3),
		Rho_yr ~ dlkjcorr(3)
		),
	data=d, 
	iter=4000, warmup=1000, chains=4, cores=4)



m2.dgs.noslope <- map2stan(
	alist(
		MidDay ~ dnorm(mu, sigma),
		mu ~ a + a_cu[cu_id] + a_yr[yr_id] + bs*fl_t + bd*date_t 
			+ bg*growth_t + ba*age_id,
		a_cu[cu_id] ~ dnorm(0, sigma_cu),
		a_yr[yr_id] ~ dnorm(0, sigma_yr),
		a ~ dnorm(15,50),
		ba ~ dnorm(15, 50),
		bs ~ dnorm(0,1),
		bg ~ dnorm(0,1),
		bd ~ dnorm(0,1),
		sigma ~ dcauchy(0,3),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3)
		),
	data=d, 
	iter=3000, warmup=1000, chains=4, cores=3)


m2.dgs.oneslope <- map2stan(
	alist(
		#likelihood
		MidDay ~ dnorm(mu, sigma),
		
		#linear models
		mu <- A + BD*date_t + BG*growth_t + BS*fl_t + BA*age_id,
		A <- a + a_cu[cu_id] + a_yr[yr_id],
		BD <- bd + bd_yr[yr_id],
		BS <- bs + bs_yr[yr_id],
		BG <- bg + bg_yr[yr_id],
		BA <- ba + ba_yr[yr_id],

		#adaptive priors
		a_cu[cu_id] ~ dnorm(0, sigma_cu),
		c(a_yr, bd_yr, bs_yr, bg_yr, ba_yr)[yr_id] ~ dmvnorm2(0, sigma_yr, Rho_yr),

		#fixed priors
		a ~ dnorm(15,50),
		c(bs, bg, bd) ~ dnorm(0,1),
		ba ~ dnorm(15, 50),
		sigma ~ dcauchy(0,3),
		sigma_cu ~ dcauchy(0,3),
		sigma_yr ~ dcauchy(0,3),
		Rho_yr ~ dlkjcorr(3)
		),
	data=d, 
	iter=4000, warmup=1000, chains=4, cores=4)


compare(m2.dgsA, m2.dgsB, m2.dgs.noslope, m2.dgs.oneslope)

precis(m2.dgsA, depth=2)
precis(m2.dgsB, depth=2)
precis(m2.dgs.noslope,depth=2)
plot(precis(m2.dgs.oneslope,depth=2))
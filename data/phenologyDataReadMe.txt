This file defines the variables included in the dataset accompanying: Freshwater, C., Trudel, M., Beacham, T.D., 
Gauthier, S., Johnson, S.C., Neville, C.-E., Juanes, F. 2018ms. Individual variation, population-specific strategies, 
and stochastic processes shape marine migration phenologies.

Note that all otolith metrics are the mean of three counts or measurements. Posterior estimates of model parameters 
were generated with fully saturated models and 9000 iterations (see manuscript for details).

Dataset: individualOtolithData.csv
captureDate		=	Julian day of capture
year			=	sampling year
latitude		= 	latitude of purse seine set
longitude		=	longitude of purse seine set
cu			=	conservation unit individual assigned to (based on microsatellite analysis described 
				in Beacham et al. 
				2004 Trans. Am. Fish. Soc.)
riverDist		=	estimated distance in km from Fraser River estuary nursery lake 
lakeLatitude		=	latitude of nursery lake
captureFL		=	individual fork length at capture (mm)
age			=	estimated age based on otolith annuli
otolithLength		=	otolith length (um)
otolithWidth		=	otolith width (um)
entryCheckRadius	=	distance from otolith core to marine entry check (um)
totalRadius		=	distance from otolith core to otolith edge (um)
incrementCount		=	number of daily increments counted between entry check and otolith edge
oceanEntryDate		=	estimate of ocean entry date (captureDate minus incrementCount)
entryFL			=	back-calculated estimate of size at ocean entry using biological intercept model (parameters
				listed in manuscript)
meanGrowth		=	average daily growth rate ((shipFL - entryFL)/incrementCount)


Dataset: ...ModPars.csv
parameter	=	parameter estimated; index references specifc conservation units (CUs) or years
mean		=	mean estimate
stDev		=	SD of estimate
lower0.89	=	lower 89% highest posterior density interval
lower0.89	=	lower 89% highest posterior density interval
effectiveN	=	estimate of effective number of samples from posterior
rHat		=	Gelman-Rubin convergence diagnostic
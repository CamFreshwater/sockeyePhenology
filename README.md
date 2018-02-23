## Individual variation, population-specific behaviours, and stochastic processes shape marine migration phenologies

This repository contains code for an analysis of juvenile sockeye salmon migration patterns using Bayesian hierarchical models and data collected from otolith microstructure. It accompanies the following paper:

Freshwater, C., M. Trudel, T.D. Beacham, S. Gauthier, S.C. Johnson, C.-E. Neville, and F. Juanes. Individual variation, population-specific behaviours, and stochastic processes shape marine migration phenologies. Submitted manuscript.

-	`modelCompare.R` uses WAIC to identify top-ranked models that included various individual characteristics as covariates
-	`fitTopModels.R` fits the top ranked model as selected by WAIC in previous script
-	`simAnalysis.R` contains a sensitivity analysis conducted to determine how sensitive the top ranked models were to the number of individual fish collected, as well as the total number of sampling years

The following R packages will be required to complete the analyses:

		install.packages(c("rethinking", "MASS", "dplyr", "gplots", "here","car"))
 	
Note that because the hierarchical models require considerable time to converge the model selection and simulation scripts may take several hours to run.

### Data
Individual otolith data are stored in `individualOtolithData.csv` and model parameters to pass directly to simulation script are in  `modPars.csv`.
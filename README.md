## BEF_quant_scale

This repository contains code to reproduce the analysis reported in the following publication:

+ coming soon (hopefully)

The aim of the project was to develop a workflow to calculate biodiversity effects on ecosystem functioning across times and places as proposed by Isbell et al. (2018, Ecology Letters) on data from natural and semi-natural systems.

The limitation to Isbell et al.'s (2018) approach is that monoculture data are required for all species in mixtures at all times and places. In addition, the method requires knowing the initial relative yields in mixtures (always 1/n species in substitutive experiments that directly manipulate species richness as a metric of biodiversity).

To solve these problems, we model monoculture yields using Bayesian generalised linear model. We then assume a range of initial relative yields in mixture by drawing initial relative yields from a Dirichlet distribution. Using the posterior distribution from the regression and the distribution of initial relative yields, we generate distributions of Isbell et al.'s (2018) biodiversity effects.

This repository includes general code to apply Isbell et al.'s (2018) partition as only limited code was provided in the regional paper (01_partition_functions). It also includes code to simulate metacommunities with full monoculture and known initial relative yields which we used to test the potential accuracy of our workflow (02_simulation). Finally, it includes code to implement our workflow to two different empirical datasets (03_empirical_analysis).

## scripts

### 01_partition_functions

In Isbell et al.'s (2018) original paper, they did not provide general code for calculating their proposed biodiversity effects. In the script *01_isbell_2018_partition.R*, we provide generalisable code that can be used for datasets with any number of species, times and places. The key function is:

+ Isbell_2018_sampler(data, RYe_post = FALSE, N = 100, alpha_par = 4, RYe)

The *data* argument is a data.frame in the formated defined by Isbell et al. (2018, Ecology Letters):

+ column 1 - **sample**: variable specifying the unique place-time combination
+ column 2 - **place**: variable specifying the place
+ column 3 - **time**: variable specifying the time-point
+ column 4 - **species**: variable specifying the species name (all times and places must have all species names present)
+ column 5 - **M**: monoculture functioning
+ column 6 - **Y**: mixture function

The *RYe* parameter is a vector of expected relative yields for the different species in the dataset. This vector must have the same length as there are species and it must sum to one.

If you do not provide RYe data and you set *RYe_post = TRUE*, the function will draw samples from the Dirichlet distribution to use as RYe values. The number of samples to draw from the Dirichlet distribution is set using *N*. The skewness of the Dirichlet distribution is set using the *alpha_par* argument.

The other two scripts: *02_isbell_2018_partition_test_data.R* and *03_isbell_2018_partition_test_script.R* are used to digitise the examples in Isbell et al. (2018) and test our functions against those examples respectively.


### 02_simulation

To test the potential accuracy of our workflow, we simulate 1000 metacommunities using Thompson et al.'s (2020, Ecology Letters) metacommunity model. Each metacommunity is simulated with slightly different assumptions about the strength of inter and intraspecific competition, dispersal rates, niche breadth etc. This model was modified slightly to simulate all monocultures for each patch in the metacommunity across times:

+ 01_mcomsimr_simulate_MC2_function.R
+ 02_analysis_test_simulate_metacommunities.R

Using these simulated metacommunities which can be generated relatively fast on a regular desktop computer, we assume that we only have monoculture data for 30% of the mixture patches. We then model the missing monoculture data using Bayesian regression implemented in Stan. This step is very computationally intensive as 1000 separate models need to be fit and sampled. We run this script on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ missing_monoculture_glm.stan
+ 03_analysis_test_model_monocultures.R

Using the posterior distributions of the simulated monocultures, we then calculate Isbell et al.'s (2018) biodiversity effects for 100 samples of initial relative yields drawn from a Dirichlet distribution. Thus, we propagate the error of the monocultures and the initial relative yields because, we use the 100 initial relative yields to calculate biodiversity for each simulated monoculture which corresponds to one sample from the posterior distribution. Therefore, if we use 1000 samples from the posterior distribution, we are left with a 1000 * 100 biodiversity effects. This script is also run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 04_analysis_test_mono_and_initial_ra_uncertainty.R

In addition, we calculate the biodiversity effects assuming known monoculture yields and unknown initial relative yields as a comparison to see how much monoculture yield uncertainty matters. This script is also run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 05_analysis_test_initial_ra_uncertainty.R

The outputs generated from the previous two steps are a very large data file. Thus, rather than processing it on our local computer, we summarise it on a script that we run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 06_analysis_test_summarise.R

The summarised outputs are then exported onto our local computer using a powershell script:

+ 07_analysis_test_export.ps1

Finally, we analyse the summarised output to test how accurate our workflow is on simulated data. This can easily done on a local computer. This script also generates figure 2.

+ 08_analysis_test_analyse.R

### Empirical analysis (03_empirical_analysis)

These scripts are used to calculate the biodiversity on the empirical data.










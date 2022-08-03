# BEF_quant_scale

Functions for quantifying the effect of biodiversity on ecosystem function across spatial and temporal scales

Here, we develop a workflow to calculate biodiversity across times and places as proposed by Isbell et al. (2018, Ecology Letters) on data from natural and semi-natural systems.

The limitation to Isbell et al.'s (2018) approach is that monoculture data is required for all species in mixtures at all times and places. In addition, it requires knowing the initial relative yields in mixtures (always 1/n species in substitutive experiments).

To solve these problems, we model monoculture yields using Bayesian regression. We then assume a range of initial relative yields in mixture by drawing initial relative yields from a Dirichlet distribution. Using the posterior distribution from the regression and the distribution of initial relative yields, we generate posterior distributions of Isbell et al.'s (2018) biodiversity effects.

This repository includes general code to apply Isbell et al.'s (2018) partition as only limited code was provided in the regional paper (01_partition_functions). It also includes code to simulate metacommunities with full monoculture and known initial relative yields which we use to test the potential accuracy of our workflow (02_simulation). Finally, it includes code to implement our workflow to different empirical datasets (03_empirical_analysis).

### Simulation (02_simulation)

To test the potential accuracy of our method, we simulate 1000 metacommunities using Thompson et al.'s (2020, Ecology Letters) metacommunity model. Each metacommunity is simulated with slightly different assumptions about the strength of inter and intraspecific competition, dispersal rates, niche breadth etc. This model was modified slightly to simulate all monocultures for each patch in the metacommunity across times:

+ 01_mcomsimr_simulate_MC2_function.R
+ 02_analysis_test_simulate_metacommunities.R

Using these simulated metacommunities which can be generated relatively fast on a regular desktop computer, we assume that we only have monoculture data for 30% of the mixture patches. We then model the missing monoculture data using Bayesian regression implemented in Stan. This step is very computationally intensive as 1000 separete models need to be fit and sampled. We run this script on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ missing_monoculture_glm.stan
+ 03_analysis_test_model_monocultures.R

Using the posterior distributions of the simulated monocultures, we then calculate Isbell et al.'s (2018) biodiversity effects for 100 samples of initial relative yields drawn from a Dirichlet distribution. Thus, we propagate the error of the monocultures and the initial relative yields because, we use the 100 initial relative yields to calculate biodiversity for each simulated monoculture which corresponds to one sample from the posterior distribution. Therefore, if we use 1000 samples from the posterior distribution, we are left with a 1000 * 100 biodiversity effects. This script is also run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 04_analysis_test_initial_ra_uncertainty.R

The output generated from the previous step is a very large data file. Thus, rather than processing it on our local computer, we summarise it on a script that we run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 05_analysis_test_summarise.R

The summarised output is then exported onto our local computer using a powershell script:

+ 06_analysis_test_export.ps1

Finally, we analyse the summarised output to test how accurate our workflow is on simulated data.

+ 07_analysis_test_analyse.R













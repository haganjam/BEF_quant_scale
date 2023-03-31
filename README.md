## BEF_quant_scale

This repository contains code to reproduce the analysis reported in the following manuscript (currently unpublished):

+ Hagan JG, Schrofner-Brunner BB, Gamfeldt L. (in prep.). Quantifying the effects of biodiversity on ecosystem functioning across space and time in natural and semi-natural ecosystems.

The aim of the project was to develop a workflow to calculate biodiversity effects on ecosystem functioning across times and places as proposed by Isbell et al. (2018, Ecology Letters) on data from natural and semi-natural systems.

The limitation to Isbell et al.'s (2018) approach is that monoculture data are required for all species in mixtures at all times and places. In addition, the method requires knowing the initial relative yields in mixtures (always 1/n species in substitutive experiments that directly manipulate species richness as a metric of biodiversity).

To solve these problems, we model monoculture yields using Bayesian generalised linear models. We then assume a range of initial relative yields in mixture by drawing initial relative yields from a Dirichlet distribution. Using the posterior distribution from the regression and the distribution of initial relative yields, we generate distributions of Isbell et al.'s (2018) biodiversity effects.

This repository includes general code to apply Isbell et al.'s (2018) partition as only limited code was provided in the regional paper (01_partition_functions). It also includes code to simulate metacommunities with full monoculture and known initial relative yields which we used to test the potential accuracy of our workflow (02_simulation). Finally, it includes code to implement our workflow to two different empirical datasets (03_empirical_analysis).

## data

Besides the simulated data which is generated using the scripts, we analysed two empirical datasets (case study 1 and case study 2 in the manuscript). 

Data and metadata for case study 1 is available on Github in the *data/case_study_1* folder. Upon publication, it will be archived on Figshare or Dryad:

+ Plymouth_data.csv
+ Plymouth_metadata.csv

Data and metadata for the second case study 2 is currently available for peer review on ResearchBox at the following link. Upon publication is will be archived on ResearchBox:

+ https://researchbox.org/843&PEER_REVIEW_passcode=GLGJFF

This unzipped research box containing all relevant raw data to reproduce the analysis needs to be downloaded and saved into the data folder on your computer. We have written a script to facilitate this process:

+ scripts/03_empirical_analysis/02_case_study_2/02_data_cleaning/01_download_data_from_researchbox.R

If you tell R where to find the unzipped researchbox, it will unzip it and save it to the data folder in this repository.

## results

The results folder stores various outputs generated throughout the analysis that may then be called later by other scripts. We have uploaded all relevant files that we generate during the analysis. Thus, to reproduce the results reported in the manuscipt, interested users do not have to run all the scripts.

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

To test the potential accuracy of our workflow, we simulate 1000 metacommunities using Thompson et al.'s (2020, Ecology Letters) metacommunity model. Each metacommunity is simulated with slightly different assumptions about the strength of inter and intraspecific competition, dispersal rates, niche breadth etc. This model was modified slightly to simulate all monocultures for each patch in the metacommunity across times.

+ 01_mcomsimr_simulate_MC2_function.R
+ 02_analysis_test_simulate_metacommunities.R

These scripts, which can be run relatively fast on a regular desktop computer, generate two files. One contains the output from the simulated metacommunities and the other contains a set of 100 samples from the Dirichlet distribution. These files are available in the results folder:

+ results/MC_sims.rds
+ results/MC_sims_start_RA.rds

Using these simulated metacommunities, we assume that we only have monoculture data for 30% of the mixture patches. We then model the missing monoculture data using Bayesian regression implemented in Stan. This step is very computationally intensive as 1000 separate models need to be fit and sampled. We run this script on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 03_analysis_test_model_monocultures.R

The stan model that the previous step calls is:

+ analysis_test_NAmonocultureGLM.stan

Using the posterior distributions of the simulated monocultures, we then calculate Isbell et al.'s (2018) biodiversity effects for 100 samples of initial relative yields drawn from a Dirichlet distribution. Thus, we propagate the error of the monocultures and the initial relative yields because, we use the 100 initial relative yields to calculate biodiversity for each simulated monoculture which corresponds to one sample from the posterior distribution. Therefore, if we use 1000 samples from the posterior distribution, we are left with a 1000 * 100 biodiversity effects. This script is also run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 04_analysis_test_mono_and_initial_ra_uncertainty.R

In addition, we calculate the biodiversity effects assuming known monoculture yields and unknown initial relative yields as a comparison to see how much monoculture yield uncertainty matters. This script is also run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 05_analysis_test_initial_ra_uncertainty.R

The outputs generated from the previous two steps are a very large data file. Thus, rather than processing it on our local computer, we summarise it on a script that we run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 06_analysis_test_summarise.R

The summarised outputs are then exported onto our local computer using a powershell script:

+ 07_analysis_test_export.ps1

The summarised outputs are provided in the *results* folder:

+ results/BEF_output.rds
+ results/BEF_output2.rds

We then analysed the summarised outputs to test how accurate our workflow is on simulated data. This can easily done on a local computer. This script also generates figure 2.

+ 08_analysis_test_analyse.R

Finally, we plot example metacommunities for the supplementary materials:

+ 09_plot_example_metacommunities.R


### 03_empirical_analysis

### 01_case_study_1

This folder contains one script that is used to calculate Isbell et al.'s (2018) biodiversity effects on data from a macroalgae removal experiment conducted in Plymouth, United Kingdom (case study 1 in the manuscript).

The folder contains one script that cleans the data, performs that analysis and makes the figures:

+ 01_plymouth_data_BEF_effects.R


### 02_case_study_2

### 01_preparation

These scripts were used to conduct a few tasks necessary before we began the experiment such as randomising the mixture and monoculture positions on the panels. They are not important for reproducing the analysis

### 02_data_cleaning

The scripts in this folder are used to unzip the data from ResearchBox, clean the different datasets: e.g. environmental data, biomass data etc. and output cleaned versions of this data into the *data/case_study_2/data_clean* folder.

The scripts are numbered 01 to 07 and should be run sequentially. All files contain a description paragraph that provides details on what the script does.

The *cleaning_functions.R* script are helper functions that are called by various cleaning scripts.

### 03_data_analysis

The scripts in this folder are used calculate Isbell et al.'s (2018) biodiversity effects on the benthic community data from the Tjarno archipelago (case study 2) and perform any other analysis reported in the manuscript.

The first script is used to test whether our temporal replicates had similar community structure using cover values from panels photographed at the same time point. The script exports Appendix 2: Figure S6.

+ 01_mixture_coverage_analysis.R

The second script analyses the environmental heterogeneity within clusters to test whether our a priori designations of clusters as heterogeneous and homogeneous based on geographic data is still valid once we use other, measured environmental variables such as temperature, salinity and water movement. The scripts exports Appendix 2: Figure S.

+ 02_environmental_variation_analysis.R

In addition, the script outputs a file to the *results* folder that contains the multivariate dispersion data:

+ results/benthic_env_dispersion.rds

The next five scripts with the prefix *03* are scripts that we used to model the missing monocultures for the five OTUs: sp A: Barn, sp B: Bryo, sp C: Asci, sp D: Hydro and sp E: Ciona. Note that these names do not directly correspond to the names in the script. This is because we renamed 'Bumpi' to 'Asci' and 'Seasq' to 'Ciona' in the manuscript for clarity.

These scripts perform model selection using PSIS. The posterior distribution of the best model and the model object are saved in the *results* folder as:

+ results/SP_X_monoculture_posterior.rds
+ results/SP_X_model_object.rds

In addition, a plot of predicted versus observed data are outputted to later combine:

+ results/SP_X_monoculture_plot.rds

The next uses the best models that were exported previously to impute the missing monocultures for each species. This script exports Appendix 1: Figures S6 and S7.

+ 03_get_monoculture_predictions.R

In addition, the script outputs two data tiles that are used to subsequently calculate the BEF effects into the *results* folder:

+ results/benthic_BEF_data.rds
+ results/benthic_start_RA.rds

Next, we used a script that takes the data with monoculture imputations and samples from the Dirichlet distribution and calculates Isbell et al.'s (2018) biodiversity effects. This script is computationally intensive and was run on a computer cluster (Albiorix: http://mtop.github.io/albiorix/).

+ 05_calculate_BEF_effects.R

We then used a powershell script to export the BEF effects to our local computer:

+ 06_benthic_analysis_export.ps1

The outputted data file was then saved into the *results* folder:

+ benthic_BEF_effects.rds

THe next script then analyses the BEF effects that we calculated and plots Figure 6 and Figure 7 from the main text:

+ 07_analyse_BEF_effects.R

We use the next script to plot the raw biomass data of mixtures and monocultures for the benthic community data: Appendix 1: Figure S8.

+ 08_plot_raw_data.R

Finally, we plot the fit of the models that we used to impute the missing monocultures: Appendix 1: Figure S5.

+ 09_plot_monoculture_model_fits.R

### renv

This project uses the renv R-package for package management. The .lock file contains all relevant information on the packages and the versions of those packages that were used in this project. To reproduce this analysis, users should install the renv R package:

install.packages("renv")
Then run the following code in the console:

renv::restore()
This will create a local copy of all relevant package versions that were used to perform these analyses.

Importantly, we were unable to add the *rethinking* package to the .lock file. Some models depend on this package. Users that would like to reproduce the results should consult the following link for download instructions (https://github.com/rmcelreath/rethinking).










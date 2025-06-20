
# BEF_quant_scale

This repository contains code to reproduce the analysis reported in the following manuscript (currently published as a preprint: https://www.researchsquare.com/article/rs-3249429/v1):

+ Hagan JG, Schrofner-Brunner B, Gamfeldt L. (in prep.). Quantifying how biodiversity affects ecosystem functioning across space and time in natural marine ecosystems.

The aim of the project was to develop a workflow to calculate biodiversity effects on ecosystem functioning across times and places as proposed by Isbell et al. (2018, Ecology Letters) on data from natural and semi-natural systems. The limitation to Isbell et al.'s (2018) approach is that monoculture data are required for all species in mixtures at all times and places. In addition, the method requires knowing the initial relative yields in mixtures (always $\frac{1}{n}$ species in substitutive experiments that directly manipulate species richness as a metric of biodiversity).

To solve these problems, we model monoculture yields using Bayesian multilevel generalised linear models. We then assume a range of initial relative yields ($RY_E$ values) in mixture by drawing initial relative yields from a Dirichlet distribution. Using the posterior distribution from the regression and the distribution of initial relative yields, we generate distributions of Isbell et al.'s (2018) biodiversity effects.

This repository includes general code to apply Isbell et al.'s (2018) partition as only limited code was provided in the regional paper (`01_partition_functions`). It also includes code to simulate metacommunities with full monoculture and known initial relative yields ($RY_E$ values) which we used to test the accuracy of our workflow (`02_simulation`). Finally, it includes code to implement our workflow to two different empirical datasets (`03_empirical_analysis`).

In addition, once the code has been run, the manuscript and supplement are automatically generated from *Quarto* documents found in the `manuscript` folder (make sure that the `make-tables.R` script is run before rendering the `manuscript.qmd` and `supplement.qmd` files).

## Running the code

To reproduce the analyses reported in the manuscript, download this repository to your computer. This can be done in two ways:

1. with git

in the Terminal:

cd path/to/local/folder

(on your computer - the folder were you want the repository to live) command on Windows might differ.

https://github.com/haganjam/BEF_quant_scale.git

This should download the directory.

2. without git

If you don't have git installed, you can download the repository as zip file and save it locally.

--> Clone or download (green button top right) --> Download Zip

then save and extract the zip where you want the directory to be saved on your computer. To run the code correctly, it is important to create a R-Project. In R-Studio go to File > New Project... > Existing Directory > Choose the extracted directory.

## data

Besides the simulated data which is generated using the scripts, we analysed two empirical datasets (case study 1 and case study 2 in the manuscript). 

Data and metadata for case study 1 are available on Github in the [**data/case_study_1**](https://github.com/haganjam/BEF_quant_scale/tree/main/data/case_study_1) folder. Upon publication, it will be archived on Figshare or Dryad:

+ Plymouth_data.csv
+ Plymouth_metadata.csv

Data and metadata for the second case study 2 are currently available for peer review on [ResearchBox](https://researchbox.org) at the following link. Upon publication it will be archived on ResearchBox:

+ https://researchbox.org/843&PEER_REVIEW_passcode=GLGJFF

This unzipped research box containing all relevant raw data to reproduce the analysis needs to be downloaded and saved into the data folder on your computer. We have written a script to facilitate this process:

+ `scripts/03_empirical_analysis/02_case_study_2/02_data_cleaning/01_download_data_from_researchbox.R`

If you tell R where to find the unzipped researchbox, it will unzip it and save it to the data folder in this repository.

## results

The results folder stores various outputs generated throughout the analysis that may then be called later by other scripts.

## scripts

### 01_partition_functions

In Isbell et al.'s (2018) original paper, they did not provide general code for calculating their proposed biodiversity effects. In the script `01_isbell_2018_partition.R`, we provide generalisable code that can be used for datasets with any number of species, times and places. The key function is:

+ `Isbell_2018_sampler(data, RYe_post = FALSE, N = 100, alpha_par = 4, RYe)`

The *data* argument is a data.frame in the format defined by Isbell et al. (2018, Ecology Letters):

+ column 1 - **sample**: variable specifying the unique place-time combination
+ column 2 - **place**: variable specifying the place
+ column 3 - **time**: variable specifying the time-point
+ column 4 - **species**: variable specifying the species name (all times and places must have all species names present)
+ column 5 - **M**: monoculture functioning
+ column 6 - **Y**: mixture function

The *RYe* parameter is a vector of expected relative yields ($RY_E$ values) for the different species in the dataset. This vector must have the same length as there are species and it must sum to one.

If you do not provide $RY_E$ data and you set *RYe_post = TRUE*, the function will draw samples from the Dirichlet distribution to use as $RY_E$ values. The number of samples to draw from the Dirichlet distribution is set using *N*. The skewness of the Dirichlet distribution is set using the *alpha_par* argument.

The other two scripts: `02_isbell_2018_partition_test_data.R` and `03_isbell_2018_partition_test_script.R` are used to digitise the examples in Isbell et al. (2018) and test our functions against those examples respectively. This is to make sure that the functions we wrote correctly calculate the effects as per Isbell et al.'s (2018) proposed method.

### 02_simulation

To test the potential accuracy of our workflow, we simulate 1000 metacommunities using Thompson et al.'s (2020, Ecology Letters) metacommunity model. Each metacommunity is simulated with slightly different assumptions about the strength of inter and intraspecific competition, dispersal rates, niche breadth etc. This model was modified slightly to simulate all monocultures for each patch in the metacommunity across times. These scripts should be run in order (i.e. 01 to 07) and can be run on a regular desktop computer in less than an hour.

### 03_empirical_analysis

#### 01_case_study_1

This folder contains one script that is used to calculate Isbell et al.'s (2018) biodiversity effects on data from a macroalgae removal experiment conducted in Plymouth, United Kingdom (case study 1 in the manuscript).

The folder contains one script that cleans the data, performs that analysis and makes the figures:

+ `01_plymouth_data_BEF_effects.R`

#### 02_case_study_2

##### 01_preparation

These scripts were used to conduct a few tasks necessary before we began the experiment such as randomising the mixture and monoculture positions on the panels. They are not important for reproducing the analysis.

##### 02_data_cleaning

The scripts in this folder are used to unzip the data from ResearchBox, clean the different datasets: e.g. environmental data, biomass data etc. and output cleaned versions of this data into the `data/case_study_2/data_clean` folder. The scripts are numbered 01 to 07 and should be run sequentially. All files contain a description paragraph that provides details on what the script does. The `cleaning_functions.R` script are helper functions that are called by various cleaning scripts.

##### 03_data_analysis

The scripts in this folder are used to calculate Isbell et al.'s (2018) biodiversity effects on the benthic community data from the Tjarno archipelago (case study 2) and perform any other analysis reported in the manuscript. As previously, the scripts need to be run in order (i.e. 01 to 08). In this folder, we also include the Stan model code that was used.

## manuscript

Once all scripts have been run, all relevant figures and table data will have been produced. The manuscript folder then contains all the necessary *Quarto* documents to reproduce the manuscript and supplement. However, the script: `make-tables.R` needs to be run before rendering the *Quarto* documents.

## renv

This project uses the renv R-package for package management. The .lock file contains all relevant information on the packages and the versions of those packages that were used in this project. To reproduce this analysis, users should install the renv R package:

install.packages("renv")
Then run the following code in the console:

renv::restore()
This will create a local copy of all relevant package versions that were used to perform these analyses.

Importantly, we were unable to add the *rethinking* package to the .lock file. Some models depend on this package. Users that would like to reproduce the results should consult the following link for download instructions (https://github.com/rmcelreath/rethinking).

### LICENCE

This file contains information about the licence under which this repository is published.







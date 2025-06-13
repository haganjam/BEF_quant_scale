#'
#' @title: Calculate the BEF effects using the imputed monoculture data
#' 
#' @description: This script executes a background job that calculates 
#' Isbell et al.'s (2018) biodiversity effects on the benthic marine fouling 
#' communities using all the variation in the modelled monocultures along with 
#' the uncertainty from the Dirichlet distribution. On my Mac M1 it took 30
#' minutes to run using 8 cores. If you do not have 8 cores, you will have to
#' modify the job code (in jobs folder)
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#'

# re-run the job: FALSE to prevent initiating the job which takes a while to run
run_calculate_BEF_effects_job <- FALSE

# run the power analysis as a background job
if (run_calculate_BEF_effects_job) {
  rstudioapi::jobRunScript(
    path = here::here("scripts/03_empirical_analysis/02_case_study_2/03_data_analysis/jobs/calculate_BEF_effects_job.R"),
    name = "calculate_BEF_effects",
    encoding = "unknown",
    workingDir = NULL,
    importEnv = FALSE,
    exportEnv = "")
}

### END

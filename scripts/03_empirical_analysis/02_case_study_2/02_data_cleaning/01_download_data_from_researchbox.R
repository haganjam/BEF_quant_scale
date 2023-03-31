#'
#' @title: Unzip the ResearchBox data and save it in the current directory
#' 
#' @description: 
#' 
#' 1. Download the data from https://researchbox.org/843&PEER_REVIEW_passcode=GLGJFF
#' 2. Run the script and select the unzipped file from ResearchBox 
#' 3. The script will unzip the fill and add it to the project directory
#' 3. You can also unzip and add it manually to "[PROJECT]/data/case_study_2"
#'
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load relevant libraries
require(here)

# unzip the file
unzip(zipfile = file.choose(),
      exdir = here("data/case_study_2"))

### END
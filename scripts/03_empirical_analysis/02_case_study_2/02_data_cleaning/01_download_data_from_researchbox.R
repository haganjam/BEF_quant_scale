#'
#' @title: Unzip the ResearchBox data and save it in the current directory
#' 
#' @description: 
#' 
#' 1. Go to the following link: https://researchbox.org/843&PEER_REVIEW_passcode=GLGJFF
#' 2. At the bottom of the page, select Download Entire Box
#' 3. Run this script and select the unzipped file that you downloaded (e.g. from downloads folder or wherever it was saved)
#' 4. The script will unzip the fill and add it to the project directory
#'
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load relevant libraries
require(here)

# unzip the file
unzip(zipfile = file.choose(),
      exdir = here("data/case_study_2"))

### END
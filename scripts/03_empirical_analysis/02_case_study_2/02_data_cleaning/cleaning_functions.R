#'
#' @title: Cleaning functions to check for directories
#' 
#' @authors: Benedikt Schrofner-Brunner
#' 

check.dirs = function(){
  
  library(here)
  
  # check that the correct folder is present
  if(! dir.exists(here("data/case_study_2/ResearchBox 843"))){
    print("download the ResearchBox contents and run cleaning script 1")
  }
  # output the cleaned csv file into the analysis data folder
  if(!dir.exists("data/case_study_2/data_clean")){ 
    dir.create("data/case_study_2/data_clean") 
  }
}

get.data.filenames = function(section.name){

  library(dplyr)
  library(here)
  
  # load the .csv file with the table of contents from ResearchBox
  table_of_contents <- read_csv(here("data/case_study_2/ResearchBox 843/Table of Contents ResearchBox 843 (peer-review).csv"))

  # get file list from researchbox dir
  table_of_contents <- 
    table_of_contents %>% 
    filter(type=="Data", section == section.name)
  
  return(table_of_contents$`file name`)
  
}



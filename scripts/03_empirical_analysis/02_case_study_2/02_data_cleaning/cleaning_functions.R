#'
#' @title: Cleaning functions to check for directories
#' 
#' @authors: Benedikt Schrofner-Brunner
#' 

check.dirs = function(){
  
  library(here)
  
  # This function checks if the folders from the Researchbox are presents, so the cleaning scripts can run.
  
  if(! dir.exists(here("data/case_study_2/ResearchBox 843"))){
    print("download the ResearchBox contents and run cleaning script 1")
  }
  # output the cleaned csv file into the analysis data folder
  if(!dir.exists("data/case_study_2/data_clean")){ 
    dir.create("data/case_study_2/data_clean") 
  }
}

#This function returns a list of filenames from the Table of Contents from Research box
get.data.filenames = function(section.name){

  library(dplyr)
  library(here)
  
  # load the .csv file with the table of contents from ResearchBox public and peer review
  
    # File names
  file_name1 <- "Table of Contents ResearchBox 843 (peer-review).csv"
  file_name2 <- "Table of Contents ResearchBox 843 (public).csv"
  file_path <- here("data/case_study_2/ResearchBox 843")
  
  if (file.exists(file.path(file_path, file_name1))) {
    
    table_of_contents <- read_csv(file.path(file_path, file_name1))
    
  } else if (file.exists(file.path(file_path, file_name2))) {
    # If the first file doesn't exist, check if the second file exists
    table_of_contents <- read_csv(file.path(file_path, file_name2))
    
  } else {stop("Table of Contents file cannot be found.")}
  
  
  # get file list from researchbox dir
  table_of_contents <- 
    table_of_contents %>% 
    filter(type=="Data", section == section.name)
  
  return(table_of_contents$`file name`)
  
}



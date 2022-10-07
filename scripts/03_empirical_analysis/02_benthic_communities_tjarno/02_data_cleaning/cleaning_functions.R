

check.dirs = function(){
  
  # check that the correct folder is present
  if(! dir.exists(here("data/benthic_communities_tjarno_data/ResearchBox 843"))){
    print("download the ResearchBox contents and run cleaning script 1")
  }
  # output the cleaned csv file into the analysis data folder
  if(!dir.exists("data/benthic_communities_tjarno_data/data_clean")){ 
    dir.create("data/benthic_communities_tjarno_data/data_clean") 
  }
}

get.data.filenames = function(section.name){

  table_of_contents <- read_csv(here("data/benthic_communities_tjarno_data/ResearchBox 843/Table of Contents ResearchBox 843 (peer-review).csv"))

    #get file list from researchbox dir
  table_of_contents = table_of_contents %>% filter(type=="Data",
                                                   section==section.name)
  table_of_contents$`file name`
}







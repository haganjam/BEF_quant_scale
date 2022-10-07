


#This script will read  all logger data and create a merged table

#This script requires:
#-Files in /data/benthic_communities_data/light_temperature_logger_data

# load relevant libraries
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(stringr)
library(lubridate)




library(readr)
library(dplyr)

#Load cleaning functions
source(here("scripts/03_empirical_analysis/02_benthic_communities_tjarno/02_data_cleaning/cleaning_functions.r"))

check.dirs()
files=get.data.filenames("light_sensor_correction")




##


# check that the correct folder is present
if(! dir.exists(here("data/benthic_communities_tjarno_data/ResearchBox 843"))){
  print("download the ResearchBox contents and run cleaning script 1")
}

table_of_contents <- read_csv(here("data/benthic_communities_tjarno_data/ResearchBox 843/Table of Contents ResearchBox 843 (peer-review).csv"))

#get file list from researchbox dir
table_of_contents = table_of_contents %>% filter(type=="Data",
                                     section=="Light_Temperature_logger_data")
files=table_of_contents$`file name`







# output the cleaned csv file into the cleaned data folder
if(!dir.exists("data/benthic_communities_tjarno_data/data_clean")){ 
  dir.create("data/benthic_communities_tjarno_data/data_clean") 
}



#replece with here
setwd("data/benthic_communities_data/light_temperature_logger_data")

files=list.files()
files[]

#ignore codebooks
files = files[!grepl("CODEBOOK",files)]

#___________________________________________________________________________
### Master table ready to use

hobos <- vector("list",length = length(files)) 
# list and how long it should be
hobos


#start loop
for(i in 1:length(files)){  #loop from i to the last number length()
  print(i)
  print(files[i])
  
  
  if(nchar(files[i])==7){
    #do this if filename is short - differences between the old and new hobo
    
    temp_table <- read_excel(
      paste("data/benthic_communities_data/light_temperature_logger_data/",
            files[i],sep=""), skip = 1)
    temp_table <- temp_table[ , 1:5]
    
    colnames(temp_table) <- c("id", "date", "time", "temp", "lux")
    temp_table[c("date","time")] <- str_split_fixed(temp_table$date," ",2)
    temp_table$site <- c(substr(files[i],1,2))
    #substr to have the first two letters of the file's name -> name of the site
    temp_table <- temp_table[,c(6,1,2,3,4,5)] #reorder the columns
    
    
  } else {
    #do this if filename is long
    
    temp_table <- read_excel(paste("data/benthic_communities_data/light_temperature_logger_data/",files[i],sep=""))
    temp_table <- temp_table[ , 1:4]
    
    colnames(temp_table) <- c("id", "date", "temp", "lux")
    temp_table[c("date","time")] <- str_split_fixed(temp_table$date," ",2)
    temp_table <- temp_table[,c(1,2,5,3,4)]
    temp_table$site <- c(substr(files[i],1,2))
    temp_table <- temp_table[,c(6,1,2,3,4,5)]
    
  }
  
  # fill up the list
  hobos[[i]] <- temp_table #write the temporary table into the master one
  
} #end loop

#View(hobos[[3]])

hobos <- bind_rows(hobos) 
#join the list into one big dataframe - do.call do the same thing

hobos %>%
  filter(lux == "Logged")
#get rid of the "Logged" text in the end of the datasheet

hobos %>%
  filter(temp == "Logged")

hobos <- hobos[complete.cases(hobos), ] 
#___________________________________________________________________________
#Filter the data when is out side of the water:
#27/07 and 28/07 and 29/07
#08/08 and 9/08 starts t1
#17/08 and 18/08 starts t2
#1/09 and 2/09 starts t3


#filter the data "inside" the date and time
hobos <- 
  hobos %>%
  mutate(date_time = ymd_hms(paste(date, time, sep = " "), tz = "Europe/Stockholm") )
#mutate is going to add new variables. Join date and time to be easier to deal with that. tz = specific time zone

View(hobos)


#Filter the time between (8h-20h) 
#define time interval with date and time

d19 <- ymd_hms("2022-07-07 00:00:00", tz = "Europe/Stockholm") #START day
d20 <- ymd_hms("2022-07-07 23:00:00", tz = "Europe/Stockholm") 
d <- ymd_hms("2022-07-27 08:00:00", tz = "Europe/Stockholm") #t0
d2 <- ymd_hms("2022-07-27 20:00:00", tz = "Europe/Stockholm") 
d3 <- ymd_hms("2022-07-28 08:00:00", tz = "Europe/Stockholm") 
d4 <- ymd_hms("2022-07-28 20:00:00", tz = "Europe/Stockholm")
d5 <- ymd_hms("2022-07-29 08:00:00", tz = "Europe/Stockholm") 
d6 <- ymd_hms("2022-07-29 20:00:00", tz = "Europe/Stockholm")
d7 <- ymd_hms("2022-08-08 08:00:00", tz = "Europe/Stockholm") #t1
d8 <- ymd_hms("2022-08-08 20:00:00", tz = "Europe/Stockholm")
d9 <- ymd_hms("2022-08-09 08:00:00", tz = "Europe/Stockholm") 
d10 <- ymd_hms("2022-08-09 20:00:00", tz = "Europe/Stockholm")
d11 <- ymd_hms("2022-08-17 08:00:00", tz = "Europe/Stockholm") #t2
d12 <- ymd_hms("2022-08-17 20:00:00", tz = "Europe/Stockholm")
d13 <- ymd_hms("2022-08-18 08:00:00", tz = "Europe/Stockholm") 
d14 <- ymd_hms("2022-08-18 20:00:00", tz = "Europe/Stockholm")
d15 <- ymd_hms("2022-09-01 08:00:00", tz = "Europe/Stockholm") #t3
d16 <- ymd_hms("2022-09-01 20:00:00", tz = "Europe/Stockholm")
d17 <- ymd_hms("2022-09-02 08:00:00", tz = "Europe/Stockholm") 
d18 <- ymd_hms("2022-09-02 20:00:00", tz = "Europe/Stockholm") #END day

d21 <- ymd_hms("2022-09-08 00:00:00", tz = "Europe/Stockholm")#akward days to delete
d22 <- ymd_hms("2022-09-08 23:00:00", tz = "Europe/Stockholm")#akward days to delete
d23 <- ymd_hms("2022-09-09 00:00:00", tz = "Europe/Stockholm")#akward days to delete
d24 <- ymd_hms("2022-09-09 23:00:00", tz = "Europe/Stockholm")#akward days to delete


#find the interval of sampling (loggers were outside of the water)
#filter the opposite

hobos <- 
  hobos %>%
  filter( !(date_time > d & date_time < d2) ) %>% 
  filter( !(date_time > d3 & date_time < d4) )%>% 
  filter( !(date_time > d5 & date_time < d6) )%>% 
  filter( !(date_time > d7 & date_time < d8) )%>% 
  filter( !(date_time > d9 & date_time < d10) )%>% 
  filter( !(date_time > d11 & date_time < d12) )%>% 
  filter( !(date_time > d13 & date_time < d14) )%>% 
  filter( !(date_time > d15 & date_time < d16) )%>% 
  filter( !(date_time > d17 & date_time < d18) )%>% 
  filter( !(date_time > d19 & date_time < d20) )%>% 
  filter( !(date_time > d21 & date_time < d22) )%>% 
  filter( !(date_time > d23 & date_time < d24) )


#___________________________________________________________________________________



# write t1, t2 and t3 in the table
hobos$time <- "" #Create the column empty

hobos$time[hobos$date <= '2022-08-16'] <- 't1' #chose the column, search for the specific date, write the "names" of that
hobos$time[(hobos$date >= '2022-08-17') & (hobos$date <= '2022-08-31')] <- 't2'
hobos$time[hobos$date >= '2022-09-01'] <- 't3'


write.csv(hobos,"data/benthic_communities_data/analysis_data/light_temperatre_data_clean.csv", row.names = FALSE)



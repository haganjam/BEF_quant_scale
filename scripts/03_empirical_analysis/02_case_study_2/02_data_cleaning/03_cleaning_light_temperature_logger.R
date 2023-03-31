#'
#' @title: Clean hobo logger data
#' 
#' @description: Cleaning light temperature loggers and merging them to one file
#' This script will also remove all times where the loggers were outside of the 
#' water. E.g. while maintaining monocultures. This script will require the raw hobologger data from Researchbox
#' 
#' @authors: Lara Martins, Benedikt Schrofner-Brunner
#' 

# load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(stringr)
library(lubridate)
library(readr)

#Load cleaning functions
source(here("scripts/03_empirical_analysis/02_case_study_2/02_data_cleaning/cleaning_functions.R"))

# check the directory
check.dirs()

# get filenames from researchbox table of contents
files <- get.data.filenames("Light_Temperature_logger_data")
print(files)

# get the file names
logger.folder <- here("data/case_study_2/ResearchBox 843/Data/")

# set-up a list to gather the logger data
hobos <- vector("list",length = length(files)) 

# read all excel files
for(i in 1:length(files)) {  #loop from i to the last number length()
  
  # at each iteration, print the i and file name
  print(i)
  print(files[i])
  
  # if file name is short i.e. only seven characters
  if(nchar(files[i])==7){
    
    # read in the logger data
    temp_table <- read_excel(paste(logger.folder, files[i], sep="/"), skip = 1)
    
    # select the first five columns
    temp_table <- temp_table[ , 1:5]
    
    # rename the columns
    colnames(temp_table) <- c("id", "date", "time", "temp_C", "lux")
    
    # split the date-time column?
    temp_table[c("date","time")] <- str_split_fixed(temp_table$date, " ", 2)
    
    # get the first two letters to have the site name
    temp_table$site <- c(substr(files[i], 1, 2))
    
    # reorder the columns
    temp_table <- temp_table[,c(6,1,2,3,4,5)]
    
    # add a column specifying the logger version
    temp_table$logger_version = "old"
    
    
    # if the file is long, then:
  } else {
    
    # read the hobo logger file
    temp_table <- read_excel(paste(logger.folder, files[i], sep="/"))
    
    # get the first four columns
    temp_table <- temp_table[ , 1:4]
    
    # renakme the four columns
    colnames(temp_table) <- c("id", "date", "temp_C", "lux")
    
    # split the date variable into a date and time variable
    temp_table[c("date","time")] <- str_split_fixed(temp_table$date," ",2)
    
    # reorder the columns
    temp_table <- temp_table[,c(1,2,5,3,4)]
    
    # get the site name by extracting the first two characters
    temp_table$site <- c(substr(files[i],1,2))
    
    # reorder the columns again
    temp_table <- temp_table[,c(6,1,2,3,4,5)]
    
    # a column specifying the logger type
    temp_table$logger_version = "new"
    
  }
  
  # fill up the list
  hobos[[i]] <- temp_table
  
} 

# join the list into one big dataframe
hobos <- bind_rows(hobos) 

# check if there is any logged text at the end
hobos %>%
  filter(lux == "Logged")
hobos %>%
  filter(temp_C == "Logged")

# remove rows with missing data
hobos <- hobos[complete.cases(hobos), ] 

# filter the data when is out side of the water (i.e. remove rows)

# dates to filter:
# 27/07 and 28/07 and 29/07
# 08/08 and 9/08 starts t1
# 17/08 and 18/08 starts t2
# 1/09 and 2/09 starts t3

# filter the data "inside" the date and time
# mutate is going to add new variables. Join date and time to be easier to deal with that. tz = specific time zone
hobos <- 
  hobos %>%
  mutate(date_time = ymd_hms(paste(date, time, sep = " "), tz = "Europe/Stockholm") )

# filter the time between (8h-20h) 
# define time interval with date and time

# dates and times to filter
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

# why are these days "awkward?"
d21 <- ymd_hms("2022-09-08 00:00:00", tz = "Europe/Stockholm") # awkward days to delete
d22 <- ymd_hms("2022-09-08 23:00:00", tz = "Europe/Stockholm") # awkward days to delete
d23 <- ymd_hms("2022-09-09 00:00:00", tz = "Europe/Stockholm") # awkward days to delete
d24 <- ymd_hms("2022-09-09 23:00:00", tz = "Europe/Stockholm") # awkward days to delete

# find the interval of sampling (loggers were outside of the water)
# filter the opposite

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

# write t1, t2 and t3 in the table
hobos$time <- "" # create the column empty

# chose the column, search for the specific date, write the "names" of that
hobos$time[hobos$date <= '2022-08-10'] <- 't1' 
hobos$time[(hobos$date >= '2022-08-10') & (hobos$date <= '2022-08-22')] <- 't2'
hobos$time[hobos$date >= '2022-08-23'] <- 't3'

# get depth of the loggers
site_data <- read_csv(here("data/case_study_2/ResearchBox 843/Data/site_data.csv"))

# select the relevant columns
site_data <- site_data[ c("site_id","panel_depth_m") ]

# rename the columns
names(site_data) <- c("site","depth")

# join the hobo data to the site data
hobos <- left_join(hobos, site_data)

# correct for new and old loggers - the old one is the standard

# load the light correction data used to predict the old logger light based on the new logger light measurements
light.correction.data <- read_csv(here("data/case_study_2/data_clean/light_sensor_correction_data.csv"))

# transform the new light data
light.correction.data$new.trans <- log(light.correction.data$new+1)

# fit a linear model with old light as the respones and new light and depth as predictors
mod.logger.pred <- lm(log(old+1) ~ new.trans*depth, data= light.correction.data)
summary(mod.logger.pred)

# transform the hobo logger light measurements (ln)
hobos$light.transformed = log(hobos$lux+1)

# prediction for log+1 transformed new logger, this needs to be transformed back wiht exp(pred)-1
hobos$lux.corrected = round(exp(predict(mod.logger.pred,
                                        newdata = data.frame(new.trans = hobos$light.transformed,
                                                             depth = hobos$depth)))-1)

#correction was also done for old loggers, those values need to be overwritten with the measured values
hobos$lux.corrected[hobos$logger_version=="old"] <- hobos$lux[hobos$logger_version=="old"]

# remove unnecessary columns
hobos <- 
  hobos %>% 
  select(-lux,-light.transformed)

# write a csv with the cleaned and corrected light and temperature data from the loggers
write_csv(hobos, here("data/case_study_2/data_clean/light_temperature_logger_clean.csv"))

### END

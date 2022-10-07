#'
#' @title: Clean hobo logger data
#' 
#' @description: 
#' Cleaning light temperature loggers and merging them to one file
#' This script will also remove all times where the loggers were outside of the 
#' water. E.g. while maintaining monocultures.
#' 
#' This script will require the raw hobologger data from Researchbox
#' 
#' @authors: Lara Martins, Benedikt Schrofner-Brunner
#' 

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
source(here("scripts/03_empirical_analysis/02_benthic_communities_tjarno/02_data_cleaning/cleaning_functions.r"))


check.dirs()# check if directory exists
#Get filenames from researchbox table of contents
files=get.data.filenames("raw data")

data.folder = here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/")



#-summarize light temp data
#-gypsum data
#-abiotic field measurements





#add combined table to data/analysis_data




################################
################################
########LARA SCRIPTS############
############START###############
#############HERE###############
################################

######LARA SCRIPT 2######

#Data preparation - Abiotic measurements

# load relevant libraries
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(stringr)
library(lubridate)

setwd("~/Universidade/EstagioSuecia/RHobo")

#Load data

abiotic <- read_excel("abiotic_measurements_data.xlsx",  col_types = c("date","text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "text", "text"))

abiotic <- abiotic[-83,] #delete weird days
abiotic <- abiotic[-87,] #delete weird days
abiotic <- abiotic[-99,] #delete weird days

colnames(abiotic) <- c("date","cluster", "site", "depth", "abiotic_temp", "salinity", "oxygen", "secchi", "cloudiness", "comment")

# # write t1, t2 and t3 in the table
# abiotic$time <- "" #Create an empty column
# 
# abiotic$time[abiotic$date==ymd('2022-07-21'] <- 't1' #chose the column, search for the specific date, write the "names" of that
# abiotic$time[abiotic$date=='2022-08-01'] <- 't2'
# abiotic$time[abiotic$date=='2022-08-31'] <- 't3'

abiotic_new = abiotic %>% 
  mutate(time = case_when(
    date == ymd("2022-07-20") ~ "t1",
    date == ymd("2022-08-01") ~ "t2",
    date == ymd("2022-08-31") ~ "t3"
  ))

abiotic_new <- abiotic_new[,c(1,2,3,11,5,6,7,8,9,10,4)]
abiotic <- abiotic_new


#save the data into csv in the folder
write.csv(abiotic,"abiotic_measurement_compilation_clean.csv", row.names = FALSE)


########LARA SCRIPT 3####


#Data preparation - Gypsum data

# load relevant libraries
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(stringr)
library(lubridate)

setwd("~/Universidade/EstagioSuecia/RHobo")

#Load data
gypsum <- read_excel("gypsum_measurement_compilation.xlsx", col_types = c("text", "text", "text", "text", "numeric", "numeric", "date", "date", "date", "date", "text", "numeric", "text"))

#_____________________________________________________________________________

#Delete NAs
gypsum <- gypsum[complete.cases(gypsum$weight_after), ]


#_____________________________________________________________________________

#Separate time_deploy and time_retrieve -> only time
gypsum[c("time_out1","time_deploy")] <- str_split_fixed(gypsum$time_deploy," ",2)

gypsum[c("time_out2","time_retrieve")] <- str_split_fixed(gypsum$time_retrieve," ",2)
gypsum <- gypsum[,-14] #delete the columns that we don't want
gypsum <- gypsum[,-14] #delete the columns that we don't want

#_____________________________________________________________________________
#Join date_deploy with date_retrieve


#filter the data "inside" the date and time
gypsum <- 
  gypsum  %>%
  mutate(date_time_deploy = ymd_hms(paste(date_deploy, time_deploy, sep = " "), tz = "Europe/Stockholm") )
#mutate is going to add new variables. Join date and time to be easier to deal with that. tz = specific time zone

gypsum <- 
  gypsum  %>%
  mutate(date_time_retrieve = ymd_hms(paste(date_retrieve, time_retrieve, sep = " "), tz = "Europe/Stockholm") )

#_____________________________________________________________________________
#Subtract the time between them

gypsum_time <- tibble(started_at = gypsum$date_time_deploy, 
                      ended_at = gypsum$date_time_retrieve) #"2020-05-27 10:03:52",

gypsum_time <- gypsum_time %>% 
  mutate(started_at = as.POSIXct(started_at,format = "%Y-%m-%d %H:%M:%S"),
         ended_at = as.POSIXct(ended_at,format = "%Y-%m-%d %H:%M:%S"))

gypsum_time$ended_at - gypsum_time$started_at


gypsum_time <-  gypsum_time %>% mutate(duration = ended_at - started_at)
colnames(gypsum_time) <- c("date_time_deploy", "date_time_retrieve", "hours_in_field")

gypsum_time<- tibble(gypsum_time) #transform the table difftime into a tibble

#_____________________________________________________________________________
#Merge gypsum_time and gypsum in the same dataset

gypsum <- merge(gypsum,gypsum_time,by= c("date_time_deploy"))
gypsum


#Clean the data
gypsum <- gypsum[,-(12)] 
gypsum <- gypsum[,c(2,3,4,5,6,7,8,9,10,1,11)] #reorder the columns
colnames(gypsum_time) <- c("date_time_deploy", "date_time_retrieve", "hours_in_field")
gypsum <- gypsum[,c(1,2,3,4,5,6,7,8,10,9,11)]

#_____________________________________________________________________________

#Apply the formula percent reduction per hour -> weight difference/weight before/time
gypsum$percent <- gypsum$weight_difference/gypsum$weight_before

#hours_in_field.y are difftime class -> numeric
gypsum$hours_in_field.y <- as.numeric(gypsum$hours_in_field.y)
gypsum$reduction_hour <- gypsum$percent/gypsum$hours_in_field.y # / by time

#_____________________________________________________________________________

# Organize the time -> t1, t2 and t3

gypsum$time <- "" #Create an empty column


gypsum$date_deploy <- as.Date(gypsum$date_deploy)

gypsum$time[gypsum$date_deploy <='2022-07-21'] <- 't1' #chose the column, search for the specific date, write the "names" of that
gypsum$time[(gypsum$date_deploy>='2022-08-02') & (gypsum$date_deploy <= '2022-08-29')] <- 't2'
gypsum$time[gypsum$date_deploy>='2022-08-29'] <- 't3'



#_____________________________________________________________________________
#save the data into csv in the folder
write.csv(gypsum,"gypsum_measurement_compilation_clean.csv", row.names = FALSE)


##############LARA SCRIPT 4#####


#Data summarize - Hobos

# load relevant libraries
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(stringr)
library(lubridate)

setwd("~/Universidade/EstagioSuecia/RHobo")

#load tables
hobos <- read_csv("hobos_measurement_compilation_clean.csv")
abiotic <- read_csv("abiotic_measurement_compilation_clean.csv") 
gypsum <- read_csv("gypsum_measurement_compilation_clean.csv")


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


#___________________________________________________________________________



### Mean and Coeficient of variation

t1_temp <- hobos %>% 
  filter(date_time<d7) %>% 
  group_by(site) %>% 
  summarise(mean_temp = mean(temp), cv_temp=sd(temp)/mean(temp)) #for the t1

#CV = σ / μ 
#cbind(t1_temp,"t1") not good ideia
t1_temp$time <- "t1" #tibble first 10 rows/columns, variable type, and not weird numbers


t2_temp <- hobos %>% 
  filter(date_time>d7 & date_time<d11) %>% 
  group_by(site) %>% 
  summarise(mean_temp = mean(temp), cv_temp=sd(temp)/mean(temp)) #for the t2
t2_temp$time <- "t2"



t3_temp <- hobos %>% 
  filter(date_time<d15) %>% 
  group_by(site) %>% 
  summarise(mean_temp = mean(temp), cv_temp=sd(temp)/mean(temp)) #for the t3
t3_temp$time <- "t3"


#Bind both tibbles
stat_table <- bind_rows(t1_temp, t2_temp, t3_temp)



# Mean of temp for all days per site

ggplot(stat_table, aes(x=site, y=mean_temp)) +
  geom_point(colour="#1F78B4") +
  labs(x="Cluster",
       y="Mean temperature (ºC)",
       title = "Mean of temperature per site registered by the loggers") +
  guides(x = guide_axis(angle = 90))
#The mean looks to the temp change


# Mean of temp and coeficient variation
ggplot(stat_table, aes(x=site, y=cv_temp)) +
  geom_point(colour="#1F78B4") +
  labs(x="Cluster",
       y="Mean temperature (ºC)",
       title = "CV per site registered by the loggers") +
  guides(x = guide_axis(angle = 90))
#The coeficient variation looks to the stability of the values





### Light Summary

hobos_sum1 <- hobos %>%
  group_by(site, time) %>%
  summarise(light = mean(lux), sd = sd(lux), cv = (sd(lux)/mean(lux))) #and join the summarize of the abiotic




#___________________________________________________________________________________

#Data Combination - Hobos + abiotics-fieldwork + gypsum

### Merge the two datasets into one - Hobos + Temperature

hobos_sum <- hobos %>%
  group_by(site, time) %>%
  summarise(mean_temp = mean(temp), sd = sd(temp), cv = (sd(temp)/mean(temp))) #and join the summarize of the abiotic


data = merge(hobos_sum,abiotic,by= c("time", "site"))
data


### Merge the two datasets into one - data + gypsum



data = merge(data,gypsum,by= c("time", "site"))
data


#save the data into csv in the folder
write.csv(abiotic,"main_data_comp_abiot_fieldwork.csv", row.names = FALSE)


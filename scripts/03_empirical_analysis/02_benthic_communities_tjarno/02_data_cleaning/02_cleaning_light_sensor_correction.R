
# create lm to correct light sensor with light sensor data

# load relevant libraries
library(readr)
library(dplyr)
library(here)
library(lubridate)
library(stringr)
library(readxl)
library(tidyr)

# load cleaning functions
source(here("scripts/03_empirical_analysis/02_benthic_communities_tjarno/02_data_cleaning/cleaning_functions.r"))

# check for the correct directories
check.dirs()

# get file names from the table of contents
files <- get.data.filenames("light_sensor_correction")

# specify the folder path that contains the logger data
logger.folder <- here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/")

# let-up an output list to store the loaded files
hobos.correction <- vector("list",length = length(files)) 

# read in the different logger files
# loop from i to the last number length()
for(i in 1:length(files)) {
  
  # if file i is a .csv file
  if(str_detect(files[i],".csv")){
    
    #read logger file from old logger
    temp.table <- read.csv(paste(logger.folder, files[i], sep="/"), skip = 1)
    
    # get the first four columns from the table
    temp.table <- temp.table[1:4]
    
    # rename the columns
    colnames(temp.table) <- c("no","datetime","temperature_C","light_lux")
    
    # convert the data.frame into a tibble()
    temp.table <- tibble(temp.table)
    
    # fix the dates
    temp.table$datetime <- mdy_hms(temp.table$datetime)
    
    # add a column for the logger version
    temp.table$logger_version = "old"
    
  } else {
    
    # read logger file from new logger (i.e. an excel file)
    temp.table = read_xlsx(paste(logger.folder, files[i], sep="/"))
    
    # rename the columns
    colnames(temp.table) = c("no","datetime","temperature_C","light_lux")
    
    # add a column for the logger version
    temp.table$logger_version = "new"
    
  }
  
  # add the depth information depth
  if(str_detect(files[i],"3m")){
    
    temp.table$depth <- 3
    
  } else {
    
    temp.table$depth <- 6
    
  }
  
  # add logger ID
  temp.table$logger_id <- str_split(files[i],"_")[[1]][6]
  
  # write the temporary table to the list
  hobos.correction[[i]] <- temp.table
  
  # remove the table
  rm(temp.table)
  
}

# bind the list into a data.frame
hobos.correction <- bind_rows(hobos.correction) 
 
# fix H4 logger - this one starts too early (how was this done?)
hobos.correction[hobos.correction$logger_id=="H4.xlsx", ]$datetime <- hobos.correction[hobos.correction$logger_id=="H4.xlsx", ]$datetime+120

# remove the first five rows
hobos.correction <- hobos.correction[-(1:5), ]

# relabel the rows
hobos.correction[hobos.correction$logger_id == "H4.xlsx",]$no <- hobos.correction[hobos.correction$logger_id == "H4.xlsx",]$no-5

# check that the loggers have the same numbers of measurements
hobos.correction %>% 
  group_by(logger_id) %>% 
  summarise(n = n())

# if there are more observations than 900, replace them with NAs
hobos.correction$no[hobos.correction$no > 900] <- NA

# remove the NAs
hobos.correction <- na.omit(hobos.correction)

# check the number of measurements again
hobos.correction %>% 
  group_by(logger_id) %>% 
  summarise(n = n())


# check for variability within the logger version

# 3m new
loggers.compare <- 
  hobos.correction %>% 
  filter(logger_version == "new", depth==3)

# check for temperature correlation
temp_C <- 
  loggers.compare %>% 
  select(no, logger_id, temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)

# calculate the correlation coefficient
cor(temp_C[,-1])

# check for variability within version
light <- 
  loggers.compare %>% 
  select(no,logger_id, light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)

# calculate the correlation coefficient
cor(light[,-1])

# rename the columns
names(light) = c("no","a","b")

# calculate the rowwise average
light_new_3 <- 
  light %>% 
  rowwise() %>% 
  mutate(mean = mean(c(a, b)))


# 6m new
loggers.compare <- 
  hobos.correction %>% 
  filter(logger_version == "new", depth==6)

# check for temperature correlation
temp_C <- 
  loggers.compare %>% 
  select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C[,-1])

# check for variability within version
light <- 
  loggers.compare %>% 
  select(no, logger_id, light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light[,-1])

# rename the light columns
names(light) <- c("no","a","b")

# calculate the rowwise means
light_new_6 <- 
  light %>% 
  rowwise() %>% 
  mutate(mean=mean(c(a,b)))


# 3m old
loggers.compare <- 
  hobos.correction %>% 
  filter(logger_version == "old", depth==3)

# check for temperature correlation
temp_C <- 
  loggers.compare %>% 
  select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C[,-1])

#check for variability within version
light <- 
  loggers.compare %>% 
  select(no,logger_id,light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light[,-1])

# rename the columns
names(light) <- c("no","a","b")

# calculate the rowwise means
light_old_3 <-  
  light %>% 
  rowwise() %>% 
  mutate(mean=mean(c(a,b)))


# 6m old
loggers.compare <- 
  hobos.correction %>% 
  filter(logger_version == "old", depth == 6)

# check for temperature correlation
temp_C <- 
  loggers.compare %>% 
  select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C[,-1])

# check for variability within version
light <- 
  loggers.compare %>% 
  select(no,logger_id,light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light[,-1])

# rename the columns
names(light) <- c("no","a","b")

# calculate the rowwise means
light_old_6 <- 
  light %>% 
  rowwise() %>% 
  mutate(mean=mean(c(a,b)))


# combine old and new loggers at different depths into a single data.frame
light.prediction.data <- 
  
  tibble(
    rbind(data.frame(old = light_old_3$mean,
                     new = light_new_3$mean,
                     depth = 3),
          data.frame(
            old = light_old_6$mean,
            new = light_new_6$mean,
            depth = 6)
          )
    )

# predict old with new
summary(lm(log(old+1)~log(new+1),data=light.prediction.data))

# plot the bivariate relationship 
plot(log(light.prediction.data$old+1),log(light.prediction.data$new+1),col=light.prediction.data$depth)

# ln-transform the old logger data
light.prediction.data$old.trans <- log(light.prediction.data$old+1)

# model the new logger light data as a function of the old logger data and depth
mod.logger.pred <- lm(log(new+1) ~ old.trans * depth, data=light.prediction.data)
summary(mod.logger.pred)

# get the predicted data from the model in original units
prediction <- 
  predict(mod.logger.pred,
          newdata = data.frame(old.trans = light.prediction.data$old.trans,
                               depth = light.prediction.data$depth))
  
# convert to original units
prediction <- exp(prediction) - 1

# correlate the predicted values based on the model with the new logger data
cor(prediction,light.prediction.data$new)

# ln() transform the new logger data
light.prediction.data$new.trans <- log(light.prediction.data$new + 1)

# fit a model that predicts old logger data using the new logger data
mod.logger.pred = lm(log(old+1) ~ new.trans * depth, data=light.prediction.data)
summary(mod.logger.pred)

# get the predictions for the old logger light data
prediction <- 
  predict(mod.logger.pred,
          newdata = data.frame(new.trans = light.prediction.data$new.trans,
                               depth = light.prediction.data$depth))
 
# convert to original units
prediction <- exp(prediction) - 1

# correlate the old data with the new data
cor(prediction,light.prediction.data$old)

# write data to create model
write_csv(light.prediction.data, here("data/benthic_communities_tjarno_data/data_clean/light_sensor_correction_data.csv"))

### END

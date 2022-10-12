#create lm to correct light sensor with light sensor data

library(readr)
library(dplyr)
library(here)
library(lubridate)

#Load cleaning functions
source(here("scripts/03_empirical_analysis/02_benthic_communities_tjarno/02_data_cleaning/cleaning_functions.r"))

check.dirs()
files=get.data.filenames("light_sensor_correction")

logger.folder = here("data/benthic_communities_tjarno_data/ResearchBox 843/Data/")


hobos.correction <- vector("list",length = length(files)) 

#read in logger files
for(i in 1:length(files)){  #loop from i to the last number length()
  
  
if(str_detect(files[i],".csv")){
  #read logger file from old logger
  temp.table = read.csv(paste(logger.folder,files[i],sep=""),skip = 1)
  
  temp.table=temp.table[1:4]
  
  colnames(temp.table) = c("no","datetime","temperature_C","light_lux")
  temp.table = tibble(temp.table)
  
  #fix date
  temp.table$datetime = mdy_hms(temp.table$datetime)
  
  temp.table$logger_version = "old"

  

  
} else {
  #read logger file from new loggger
  
  temp.table = read_xlsx(paste(logger.folder,files[i],sep=""))
  colnames(temp.table) = c("no","datetime","temperature_C","light_lux")
  temp.table$logger_version = "new"


}
  #add depth
  if(str_detect(files[i],"3m")){
    temp.table$depth=3
  } else {
    temp.table$depth=6
  }
  
  #add logger ID
  temp.table$logger_id=str_split(files[i],"_")[[1]][6]
  
  hobos.correction[[i]] = temp.table
  
  rm(temp.table)
  
}

hobos.correction <- bind_rows(hobos.correction) 

#fix H4 logger - this one starts too early
hobos.correction[hobos.correction$logger_id=="H4.xlsx",]$datetime=hobos.correction[hobos.correction$logger_id=="H4.xlsx",]$datetime+120
hobos.correction = hobos.correction[-(1:5),]
hobos.correction[hobos.correction$logger_id=="H4.xlsx",]$no=hobos.correction[hobos.correction$logger_id=="H4.xlsx",]$no-5

#have same number of measurements
hobos.correction %>% group_by(logger_id) %>% summarise(n())
hobos.correction$no[hobos.correction$no>900]=NA
hobos.correction = na.omit(hobos.correction)
hobos.correction %>% group_by(logger_id) %>% summarise(n())


######check for variability within version#####
#3m new
loggers.compare = hobos.correction %>% filter(logger_version=="new",depth==3)

#Check for temperature correlation
temp_C = loggers.compare %>% select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C) # .981

#check for variability within version
light = loggers.compare %>% select(no,logger_id,light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light) # .918

names(light)=c("no","a","b")
light_new_3 = light %>% rowwise() %>% mutate(mean=mean(c(a,b)))


#6m new
loggers.compare = hobos.correction %>% filter(logger_version=="new",depth==6)

#Check for temperature correlation
temp_C = loggers.compare %>% select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C) # .996

#check for variability within version
light = loggers.compare %>% select(no,logger_id,light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light) # .866

names(light)=c("no","a","b")
light_new_6 = light %>% rowwise() %>% mutate(mean=mean(c(a,b)))


#3m old
loggers.compare = hobos.correction %>% filter(logger_version=="old",depth==3)

#Check for temperature correlation
temp_C = loggers.compare %>% select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C) # .996

#check for variability within version
light = loggers.compare %>% select(no,logger_id,light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light) # .479

names(light)=c("no","a","b")
light_old_3 = light %>% rowwise() %>% mutate(mean=mean(c(a,b)))


#6m old
loggers.compare = hobos.correction %>% filter(logger_version=="old",depth==6)

#Check for temperature correlation
temp_C = loggers.compare %>% select(no,logger_id,temperature_C) %>%
  pivot_wider(names_from = logger_id, values_from = temperature_C)
cor(temp_C) # .991

#check for variability within version
light = loggers.compare %>% select(no,logger_id,light_lux) %>%
  pivot_wider(names_from = logger_id, values_from = light_lux)
cor(light) # .944

names(light)=c("no","a","b")
light_old_6 = light %>% rowwise() %>% mutate(mean=mean(c(a,b)))


#combine old and new loggers
light.prediction.data = tibble(rbind(
  data.frame(
    old=light_old_3$mean,
    new=light_new_3$mean,
    depth=3),
  data.frame(
    old=light_old_6$mean,
    new=light_new_6$mean,
    depth=6)))

#Predict old with new
summary(lm(log(old+1)~log(new+1),data=light.prediction.data))


plot(log(light.prediction.data$old+1),log(light.prediction.data$new+1),col=light.prediction.data$depth)

#transform
light.prediction.data$old.trans = log(light.prediction.data$old+1)

mod.logger.pred = lm(log(new+1) ~ old.trans * depth, data=light.prediction.data)
summary(mod.logger.pred)

prediction = exp(
  predict(mod.logger.pred,
            newdata = data.frame(old.trans = light.prediction.data$old.trans,
                                 depth = light.prediction.data$depth)))-1

cor(prediction,light.prediction.data$new)

light.prediction.data$new.trans = log(light.prediction.data$new+1)

mod.logger.pred = lm(log(old+1) ~ new.trans * depth, data=light.prediction.data)
summary(mod.logger.pred)

prediction = exp(
  predict(mod.logger.pred,
          newdata = data.frame(new.trans = light.prediction.data$new.trans,
                               depth = light.prediction.data$depth)))-1

cor(prediction,light.prediction.data$old)

#it works better to predict old logger with new



#write data to create model
write.csv(light.prediction.data,here("data/benthic_communities_tjarno_data/data_clean/light_sensor_correction_data.csv"), row.names = FALSE)

rm(list = ls())

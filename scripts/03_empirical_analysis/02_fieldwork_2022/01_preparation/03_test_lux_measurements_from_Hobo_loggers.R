
# Examine the difference between the two-types of light loggers

# load relevant libraries
library(dplyr)
library(readr)
library(ggplot2)
library(here)
library(lubridate)
library(tidyr)

# useful regex examples
# sub("\\s+.*", "",  x$`Data Hora, GMT+02:00`[1:5])
# sub(".*\\s[0]", "",  x$`Data Hora, GMT+02:00`[1:5])
# y <- sub("\\s[P][M]", "", x$`Data Hora, GMT+02:00`[1:5])
# z <- sub("/", "-", x$`Data Hora, GMT+02:00`[1:5])
# u <- sub("/", "-", z)

# load the data from the new hobo loggers
files <- list.files(here("data"))

# subset the old files
old_files <- files[grepl(pattern = "_old.csv", files)]
old_files_depth <- gsub(pattern = "\\_.*", "", old_files)

old_files_out <- vector("list")
for(i in 1:length(old_files)) {
  
  x <- read_csv(file = paste(here("data"), "/", old_files[i], sep = ""),
                skip = 1)
  x <- x[, c(2, 3, 4)]
  names(x) <- c("date_time_CEST", "temp_C", "Lux")
  
  # parse the date-time column and convert to a datetime object
  u <- parse_date_time(x$date_time_CEST,  
                       "mdy HMS Op", tz = "Europe/Stockholm")
  u <- lubridate::as_datetime(u, tz = "Europe/Stockholm")
  
  # replace the date_time_CEST column
  x$date_time_CEST <- u
  
  # add a depth column
  x$depth <- as.numeric(old_files_depth[i])
  
  # add a old_new variable column
  x$old_new <- "old"
  
  # add a number column i.e. logger ID
  x$id <- i
  
  old_files_out[[i]] <- x
  
}


# subset the new files
new_files <- files[grepl(pattern = "_new.csv", files)]
new_files_depth <- gsub(pattern = "\\_.*", "", new_files)

new_files_out <- vector("list")
for(i in 1:length(new_files)) {
  
  x <- read_csv(file = paste(here("data"), "/", new_files[i], sep = ""))
  x <- x[, c(2, 3, 4)]
  names(x) <- c("date_time_CEST", "temp_C", "Lux")
  
  # parse the date-time column and convert to a datetime object
  u <- parse_date_time(x$date_time_CEST,  
                       "mdy HMS", tz = "Europe/Stockholm")
  u <- lubridate::as_datetime(u, tz = "Europe/Stockholm")
  
  # replace the date_time_CEST column
  x$date_time_CEST <- u
  
  # add a depth column
  x$depth <- as.numeric(new_files_depth[i])
  
  # add a old_new variable column
  x$old_new <- "new"
  
  # add a number column i.e. logger ID
  x$id <- i
  
  new_files_out[[i]] <- x
  
}

# join the old and new data together
m.dat <- 
  bind_rows(bind_rows(old_files_out),
            bind_rows(new_files_out))

# reorganise the columns
m.dat <- 
  m.dat %>%
  select(old_new, id, depth, date_time_CEST, temp_C, Lux) %>%
  arrange(old_new, id, depth, date_time_CEST)
head(m.dat)

# only consider data after 18h00
m.dat <- 
  m.dat %>%
  filter(date_time_CEST > as.POSIXct(x = "2022-06-30 18:00:00"),
         date_time_CEST < as.POSIXct(x = "2022-07-04 13:17:00")) %>%
  filter(!is.na(Lux))

m.dat %>%
  group_by(old_new, depth) %>%
  summarise(n = n())

m.dat %>%
  filter(old_new == "new", depth == 155) %>%
  View()

# remove depth 155 because there is no old comparison
m.dat <- 
  m.dat %>%
  filter(depth != 155)

# create the data.frame needed to test correlation between different light values
old.dat <- m.dat[m.dat$old_new == "old", ]
new.dat <- m.dat[m.dat$old_new == "new", ]

unique(old.dat$id)
unique(new.dat$id)

View(old.dat)
View(new.dat)

# change the relevant light names
old.dat <- 
  old.dat %>%
  rename(Lux_old = Lux)

new.dat <- 
  new.dat %>%
  rename(Lux_new = Lux)

l.dat <- 
  full_join(new.dat %>%
            select(depth, date_time_CEST, Lux_new), 
          old.dat %>%
            select(depth, date_time_CEST, Lux_old), 
          by = c("depth", "date_time_CEST"))

# plot the relationship
ggplot(data = l.dat %>% mutate(depth = as.character(depth)),
       mapping = aes(x = Lux_new, y = Lux_old, colour = depth)) +
  geom_point() +
  facet_wrap(~depth, scales = "free") +
  theme_bw()

ggplot(data = l.dat %>% mutate(depth = as.character(depth)),
       mapping = aes(x = log(Lux_new), y = log(Lux_old), colour = depth)) +
  geom_point() +
  facet_wrap(~depth) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()







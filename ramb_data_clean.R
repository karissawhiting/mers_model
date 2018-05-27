
library(knitr)
library(tidyverse)
library(rprojroot)
library(lubridate)
library(ggmap)
library(maps)
library(scales)
library(ggmap)
library(dygraphs)
library(forcats)
library(reshape2)
#library(gganimate)
#library(animation)

# Clean Data ----------------------
setwd("~/Repositories/practicum")
#mers_c2 <- read.csv("./mers_cases_rambaut_copy.csv", na.strings = c("", "NA"))
mers_c <- read.csv("./mers_cases_rambaut.csv", na.strings = c("", "NA"))


#Clean Mers Cases Data
mers_c<- dplyr::select(mers_c, country, city, exposure, onset, reported) 
mers_c<- mutate(mers_c, onset = ymd(onset), reported = ymd(reported), exposure = ymd(exposure)) 

#mers_c <- mers_c %>%
#  mutate(week_year = as.Date(cut(onset, breaks = "week", start.on.monday = FALSE))) %>%
#  mutate(mon_year =  as.Date(cut(onset, breaks = "month")))

# Subset SK data 
mers_sa <- mers_c %>%
  filter(country == "South Korea")

purrr::map(mers_sa, ~sum(is.na(.x)))

# get prior reporting rates
miss_rate<- purrr::map_df(mers_sa, ~sum(is.na(.x))/186)
purrr::map(miss_rate, ~1-.x)


purrr::map_df(mers_sa, ~sum(is.na(.x)))


# Get in time series form ----------------------
mers <- mers_sa %>% 
  mutate(onset_date = as.Date(onset, format="%m/%d/%y"),
         conf_date = as.Date(reported, format="%m/%d/%y"), 
         exp_date = as.Date(exposure, format="%m/%d/%y"))

miss_rate<- purrr::map_df(mers, ~sum(is.na(.x))/186)
purrr::map(miss_rate, ~1-.x)

mers_times <- mers %>%
  dplyr::select(exp_date, onset_date, conf_date) %>%
  mutate(exp_date = exp_date - min(onset_date, na.rm = TRUE) +1 ,
         conf_date = conf_date - min(onset_date, na.rm = TRUE) +1 ,
         onset_date = onset_date - min(onset_date, na.rm = TRUE) +1) 

save(mers_times, file = "../mers_model/data-raw/mers_times_ram.RData")

miss_rate<- purrr::map_df(mers_times, ~sum(is.na(.x))/186)
purrr::map(miss_rate, ~1-.x)

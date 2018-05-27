library(tidyverse)
library(MASS)

# Read in data
setwd("~/Repositories/mers_model")
mers <- read_csv("./data-raw/mers_cases_rambaut.csv")

sk <- mers %>% 
  filter(country == "South Korea") %>%
  dplyr::select(code, province, city, district, exposure, onset, hospitalized,
         reported, death, discharged, secondary) %>%
  mutate(onset = as.Date(onset, format="%m/%d/%y"),
         reported = as.Date(reported, format="%m/%d/%y"), 
         no_days = reported - onset)

####################
# Impute Data
hist(unclass(sk$no_days), breaks = 40)

vec<- na.omit(unclass(sk$no_days))
vec<- vec[vec != 0]
den <- density(vec)
dat <- data.frame(x = den$x, y = den$y)

fit.params  <- fitdistr(vec, dgamma, list(shape = 1, rate = 0.1), lower = .001)
x <- rgamma(100000, shape = gamma$estimate[1], rate = gamma$estimate[2])
hist(vec)

####################
# SEIC Model





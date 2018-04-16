library(tidyverse)
library(MASS)
library(fitdistrplus)

# Clean Data ----------------------
# Read in data
setwd("~/Repositories/mers_model")
mers_raw <- read_csv("./data-raw/KCDC_mers_cl.csv")

# clean data and turn to dates
mers <- mers_raw %>% 
  mutate(onset_date = as.Date(onset_date, format="%m/%d/%y"),
         conf_date = as.Date(conf_date, format="%m/%d/%y"), 
         f_exp_date = as.Date(f_exp_date, format="%m/%d/%y"),
         l_exp_date = as.Date(l_exp_date, format="%m/%d/%y"),
         dis_death_date = as.Date(dis_death_date, format="%m/%d/%y"),
         delay_on_conf =  conf_date - onset_date)

# check for na's
purrr::map(mers, ~sum(is.na(.x)))


# Onset Dates ----------------------
# impute missing onset dates for 7 cases

delay_on_conf <- mers$conf_date - mers$onset_date
sum(is.na(delay_on_conf))

delay<- as.numeric(na.omit(delay_on_conf))
hist(delay, breaks = 20)

delay[delay == 0] <- .00001 

hist(delay, breaks = 25)

fit.gamma <- fitdist(delay, distr = "gamma", 
                     method = "mme")
summary(fit.gamma)

plot(fit.gamma, breaks = 20)



# replace NA values with imputation

set.seed(7)
mers$onset_date = replace(mers$onset_date, which(is.na(mers$onset_date)), 
                        mers$conf_date[which(is.na(mers$onset_date))] - round(rgamma(7, shape = fit.gamma$estimate['shape'], rate =  fit.gamma$estimate['rate'])))


# Exposure Dates ----------------------

exp_win <- mers$l_exp_date - mers$f_exp_date

mers2<- mers %>%
  mutate(exp_win_days = l_exp_date - f_exp_date, 
         exp_date = NA)

for (i in 1:nrow(mers2)) {
  mers2$exp_date[i] = round(runif(1, min = 0, max = mers2$exp_win_days[i]))
}

mers2$exp_date2 = mers2$f_exp_date + mers2$exp_date


mers2 <- mers2 %>%
  mutate(inc = onset_date - exp_date2)

hist(as.numeric(mers2$inc), breaks = 40)


fit.gamma <- fitdist(as.numeric(na.omit(mers2$inc)), distr = "gamma", 
                     method = "mme")
summary(fit.gamma)

plot(fit.gamma, breaks = 20)

# replace NA values with inc imputation
set.seed(7)
mers2$exp_date2 = replace(mers2$exp_date2, which(is.na(mers2$exp_date2)), 
                          mers2$onset_date[which(is.na(mers2$exp_date2))] - 
                           round(rgamma(sum(is.na(mers2$exp_date2)), fit.gamma$estimate['shape'], rate =  fit.gamma$estimate['rate'])))


# turn into time series
mers_times <- mers2 %>%
  dplyr::select(id, exp_date2, onset_date, conf_date) %>%
  mutate(exp_date = exp_date2 - min(onset_date) +1 ,
         conf_date = conf_date - min(onset_date) +1 ,
         onset_date = onset_date - min(onset_date) +1) %>%
  dplyr::select(-exp_date2)



save(mers_times, file = "./data-raw/mers_times.RData")


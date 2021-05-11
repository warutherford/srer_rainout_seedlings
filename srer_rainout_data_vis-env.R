# USDA Rainout - Santa Rita Experimental Range
# Data visualizations-site environmental data
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-05-10


# load packages
# Load packages
library(tidyverse)
library(lubridate)
library(forcats)
library(vroom)
library(ggpubr)
library(glmmTMB)
library(dplyr)

### soil temp

stemplight <- vroom("Data/site-env-data/soil-temp/all_soil_temp_light.csv",
               col_select = -c(1),
               col_types = c(.default = "d",
                             datetime = "T",
                             precip = "f",
                             location = "f",
                             section = "f",
                             clip = "f"))
str(stemplight)


# pull out only 5cm soil temp data
stemp <- stemplight %>% 
  dplyr::select(-lux) %>% 
  filter(location == 'BelowTemp') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  group_by(precip, location, clip) %>% 
  summarise(tempC_mean = mean(soiltempC),
            tempC_se = sd(soiltempC)/sqrt(n()))

# pull out only surface air temp
atemp <- stemplight %>% 
  dplyr::select(-lux) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  group_by(precip, location, clip) %>% 
  summarise(atempC_mean = mean(soiltempC),
            atempC_se = sd(soiltempC)/sqrt(n()))

# time series of 5cm soil temp
stemp_series <- stemplight %>% 
  dplyr::select(-lux) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  filter(location == 'BelowTemp') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  mutate(hour = hour(datetime),
         week = week(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>%
  summarise(dailyavg = mean(soiltempC))

stemp_series %>% 
  group_by(precip, clip, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg, color = clip, linetype = year))+
  geom_line(size = 2)+
  xlim(0,52)+
  labs(y = "Daily 5-cm Soil Temperature (°C)",
       x = "Week",
       color = "Grazing/Clipping",
       linetype = "Year") +
  labs_pubr() +
  facet_wrap(~precip+year) +
  theme_pubr(legend = "bottom")

# time series of air surface soil temp
atemp_series <- stemplight %>% 
  dplyr::select(-lux) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00'))

# make bad values, NAs
atemp_fix <- atemp_series %>% 
  filter(precip == 'Ambient' & section == 'E1') %>% 
  filter((datetime >= ymd('2019-05-06') & datetime <= ymd('2019-08-01'))) %>% 
  mutate(soiltempC = replace(soiltempC, soiltempC>0, NA))

atemp_comp <- full_join(atemp_series, atemp_fix) %>% 
  mutate(hour = hour(datetime),
         week = week(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>% 
  summarise(dailyavg = mean(soiltempC))

# plot time series with missing data in ambient clipped from 2019
surface_temp_fig <- atemp_comp %>%
  group_by(precip, clip, week, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg, color = clip, linetype = year)) +
  geom_line(size = 2) +
  scale_x_continuous(breaks=seq(0, 52, 4)) +
  #xlim(NA, 52) +
  labs(y = "Daily Soil Suface Air Temperature (°C)",
       x = "Week",
       color = "Grazing/Clipping",
       linetype = "Year") +
  labs_pubr() +
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom")

surface_temp_fig


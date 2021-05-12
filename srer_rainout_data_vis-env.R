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
library(hydroTSM)
library(dplyr)

### Station PPT and Air Temp
ppt_temp <- vroom("Data/site-env-data/ppt-temp/all_air_ppt_temp.csv",
                    col_types = c(.default = "d",
                                  datetime = "T",
                                  pptmm = "n",
                                  tempC = "n",
                                  ids = "f"))
str(ppt_temp)


ppt_temp_series <- ppt_temp %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  extract("ids", into = c("site"), "(.{4})", remove = TRUE) %>% 
  mutate(site = as.factor(site),
         pptmm = replace(pptmm, is.na(pptmm), 0))

glimpse(ppt_temp_series)
summary(ppt_temp_series)

test <- ppt_temp_series %>%
  distinct(datetime, .keep_all = TRUE) %>% 
  mutate(day = as.factor(yday(datetime)),
         month = as.factor(month(datetime)),
         year = as.factor(year(datetime))) %>%
  group_by(site, day, year) %>% 
  summarise(ppt_measure = max(pptmm)) %>% 
  arrange(site , year, day)

test_diff <- test %>%
  filter(ppt_measure != 0) %>% 
  group_by(site) %>% 
  mutate(ppt_diff = (ppt_measure - dplyr::lag(ppt_measure, default = last(ppt_measure)))) %>% 
  mutate(ppt_diff = replace(ppt_diff, ppt_diff < 0, 0))

test_zero <- test %>% 
  filter(ppt_measure == 0)

test_comp <- rbind(test_diff, test_zero)
  
  
test_comp_all <- test_comp %>% mutate(ppt_diff = replace(ppt_diff, is.na(ppt_diff), 0)) %>% 
  arrange(site, year, day)

summary(test_diff)

test_avg <- test_comp_all %>% 
  group_by(day, year) %>% 
  summarise(ppt_mean = mean(ppt_diff))
  

summary(test_avg)
glimpse(test_avg)

test_avg %>% 
  filter(year != '2020') %>% 
  #group_by(day, week, month, year) %>% 
ggplot(aes(x = day, y = ppt_mean, color = year)) +
  geom_bar(stat = "identity") +
  facet_wrap(~year, ncol = 3)

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

# make sure date range is within seedling observation date range
stemp_stats <- stemplight %>% 
  dplyr::select(-lux) %>% 
  filter(location == 'BelowTemp') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE)


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

# test sig diff of daily below ground soil temp between treatments
# data normal?
hist(stemp_series$dailyavg) # loods good

below_anova <- aov(dailyavg ~ precip*clip*year, data = stemp_series)

summary(below_anova)

post.hoc.below <- emmeans::emmeans(below_anova, specs = ~precip+year)
post.hoc.below

# get lettering report on post-hoc test
post.hoc.letters.below <- cld(post.hoc.below, Letters = letters, covar = T)
post.hoc.letters.below

# figure of 5-cm soil temp time series
below_soil_temp_fig <- stemp_series %>% 
  group_by(precip, clip, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg, color = clip, linetype = year))+
  geom_line(size = 2)+
  #geom_smooth() +
  scale_x_continuous(breaks=seq(0, 52, 4)) +
  #xlim(0,52)+
  labs(y = "Daily 5-cm Soil Temperature (°C)",
       x = "Week",
       color = "Grazing/Clipping",
       linetype = "Year") +
  labs_pubr() +
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom")

below_soil_temp_fig

ggsave(filename = "Figures_Tables/5cm_soil_temp.tiff",
       plot = below_soil_temp_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Air Soil Surface temp

# time series of air surface soil temp, make sure time is within seedling observations
atemp_series <- stemplight %>% 
  dplyr::select(-lux) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00'))

# Corrupted sensor data make NA
atemp_fix <- atemp_series %>% 
  filter(precip == 'Ambient' & section == 'E1') %>% 
  filter((datetime >= ymd('2019-05-06') & datetime <= ymd('2019-08-01'))) %>% 
  mutate(soiltempC = replace(soiltempC, soiltempC>0, NA))

# join full data set with NAs, extract daily averages for graphing
atemp_comp <- full_join(atemp_series, atemp_fix) %>% 
  mutate(hour = hour(datetime),
         week = week(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>% 
  summarise(dailyavg = mean(soiltempC))

# test sig diff of daily surface soil temp between treatments
# data normal?
hist(atemp_comp$dailyavg) # normal

atemp_aov <- aov(dailyavg ~ precip*clip*year, data = atemp_comp)

summary(atemp_aov)

post.hoc.surf <- emmeans::emmeans(atemp_aov, specs = ~precip+year)
post.hoc.surf

# get lettering report on post-hoc test
post.hoc.letters.surf <- cld(post.hoc.surf, Letters = letters, covar = T)
post.hoc.letters.surf

# plot time series with missing data in ambient clipped from 2019
surface_temp_fig <- atemp_comp %>%
  group_by(precip, clip, week, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg, color = clip, linetype = year)) +
  geom_line(size = 2) +
  #geom_smooth() +
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

ggsave(filename = "Figures_Tables/surface_soil_temp.tiff",
       plot = surface_temp_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Soil Surface Light
# extract light data
light_stats <- stemplight %>% 
  drop_na(lux) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00') & 
           lux > 0) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  group_by(precip, location, clip) %>% 
  summarise(light_mean = mean(lux),
            light_se = sd(lux)/sqrt(n()))

# make sure date range is within seedling observation date range
light <- stemplight %>% 
  dplyr::select(-soiltempC) %>% 
  drop_na(lux) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00') & 
           lux > 0) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  group_by(precip, location, clip)

# Corrupted sensor data make NA
light_fix <- light %>% 
  filter(precip == 'Ambient' & section == 'E1') %>% 
  filter((datetime >= ymd('2019-05-06') & datetime <= ymd('2019-08-01'))) %>% 
  mutate(lux = replace(lux, lux > 0, NA))

# join full data set with NAs, extract daily averages for graphing
light_comp <- full_join(light, light_fix) %>% 
  mutate(hour = hour(datetime),
         week = week(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>% 
  summarise(dailyavg_light = mean(lux))

# test sig diff of daily below ground soil temp between treatments
# data normal?
hist(light_series$dailyavg_light) # loods good

light_anova <- aov(dailyavg_light ~ precip*clip*year, data = light_series)

summary(light_anova)

post.hoc.light <- emmeans::emmeans(light_anova, specs = ~precip+clip)
post.hoc.light

# get lettering report on post-hoc test
post.hoc.letters.light <- cld(post.hoc.light, Letters = letters, covar = T)
post.hoc.letters.light

# figure of light time series
light_fig <- light_comp %>% 
  group_by(precip, clip, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg_light, color = clip, linetype = year))+
  geom_point(size = 2)+
  scale_x_continuous(breaks=seq(0, 52, 4)) +
  scale_y_continuous(labels = scales::scientific) + 
  geom_smooth() +
  #xlim(0,52)+
  labs(y = "Daily Light (lux)",
       x = "Week",
       color = "Grazing/Clipping",
       linetype = "Year") +
  labs_pubr() +
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom")

light_fig

ggsave(filename = "Figures_Tables/surface_light.tiff",
       plot = light_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

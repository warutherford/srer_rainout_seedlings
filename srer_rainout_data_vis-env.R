# USDA Rainout - Santa Rita Experimental Range
# Data visualizations-site environmental data
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-05-10

# Load packages
library(tidyverse)
library(lubridate)
library(forcats)
library(vroom)
library(ggpubr)
library(glmmTMB)
library(dplyr)
library(multcomp)
library(psych)

### Santa Rita Desert Grassland Station Rain Can/DESGR
ppt_srer <- vroom("Data/site-env-data/precip_srer_DESGR.csv",
                  col_types = c(.default = "n",
                                STATION = "c"))

glimpse(ppt_srer)

ppt_srer_avg <- ppt_srer %>%
  group_by(YEAR) %>% 
  rowwise() %>% 
  summarise(ppt_mean = mean(c_across(JAN:DEC), na.rm = TRUE)) %>% 
  mutate(ppt_yr_in = ppt_mean/100) %>% 
  mutate(ppt_yr_mm = ppt_yr_in*25.4) %>% 
  ungroup()

summary(ppt_srer_avg)


### Rainout Station PPT and Air Temp
ppt_temp <- vroom("Data/site-env-data/ppt-temp/all_air_ppt_temp.csv",
                    col_types = c(.default = "d",
                                  datetime = "T",
                                  pptmm = "n",
                                  tempC = "n",
                                  ids = "f"))
str(ppt_temp)


ppt_temp_series <- ppt_temp %>% 
  filter(datetime >= ymd_hms('2017-01-01 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  tidyr::extract("ids", into = c("site"), "(.{4})", remove = TRUE) %>% 
  mutate(site = as.factor(site),
         pptmm = replace(pptmm, is.na(pptmm), 0))

glimpse(ppt_temp_series)
summary(ppt_temp_series)

ppt_parse <- ppt_temp_series %>%
  distinct(datetime, .keep_all = TRUE) %>% 
  mutate(day = as.factor(yday(datetime)),
         week = week(datetime),
         month = as.factor(month(datetime, label = TRUE)),
         year = as.factor(year(datetime)),
         monthday = as.factor(mday(datetime))) %>%
  filter(year != 2020) %>% 
  group_by(site, day, week, monthday, month, year) %>% 
  summarise(ppt_measure = max(pptmm)) %>% 
  arrange(site , year, month, week, day)

ppt_lag <- ppt_parse %>%
  filter(ppt_measure != 0) %>% 
  group_by(site) %>% 
  mutate(ppt_diff = (ppt_measure - dplyr::lag(ppt_measure, default = 0))) %>% 
  mutate(ppt_diff = replace(ppt_diff, ppt_diff < 0, NA))
  
ppt_filt <- ppt_lag %>%   
  filter(is.na(ppt_diff)) %>% 
  mutate(ppt_diff = ppt_measure)

ppt_zero <- ppt_parse %>% 
  filter(ppt_measure == 0)

ppt_comp <- rbind(ppt_lag, ppt_zero, ppt_filt)
  
  
ppt_comp_all <- ppt_comp %>% 
  mutate(ppt_diff = replace(ppt_diff, is.na(ppt_diff), 0)) %>% 
  arrange(site, year, month, week, day)

ppt_summary <- ppt_comp_all %>% 
  group_by(day, week, monthday, month, year) %>% 
  summarise(ppt = mean(ppt_diff)) %>% 
  mutate(ppt_in = ppt/100,
         ppt_mm = ppt*0.254) %>% 
  rowid_to_column('id') %>% 
  filter(id != 1) %>% 
  dplyr::select(-id)

ppt_summary_site <- ppt_comp_all %>% 
  group_by(day, week, monthday, month, year, site) %>% 
  summarise(ppt = mean(ppt_diff)) %>% 
  mutate(ppt_in = ppt/100,
         ppt_mm = ppt*0.254) %>% 
  rowid_to_column('id') %>% 
  filter(id != 1) %>% 
  dplyr::select(-id) 

ppt_start <-  ppt_comp_all %>% 
  group_by(day, week, monthday, month, year) %>% 
  summarise(ppt = mean(ppt_diff)) %>% 
  mutate(ppt_in = ppt/100,
         ppt_mm = ppt*0.254) %>%
  rowid_to_column("id") %>% 
  filter(id == 1) %>% 
  mutate(ppt = 0,
         ppt_in = 0,
         ppt_mm = 0) %>% 
  dplyr::select(-id)

ppt_final <- full_join(ppt_summary, ppt_start) %>% 
  arrange(year, month, week, day)

ppt_final %>% 
  group_by(year) %>% 
  summarise(ppt_sum = sum(ppt_mm))

# PPT figure all years
ppt_all_fig <- ppt_final %>% 
  group_by(week, year) %>% 
  ggplot(aes(x = week, y = ppt_mm, fill = year)) +
  geom_col() +
  #geom_hline(yintercept = 34.25, color = "darkgreen", size = 1.5) +
  scale_fill_manual(values = c("#b0e8f5", "#169cf0", "#0033FF")) +
  scale_x_continuous(breaks=seq(0, 52, 4)) +
  scale_y_continuous(breaks = seq(0, 110, 10), expand = c(0.01,0)) +
  #xlim(0,52)+
  labs(y = "Precipitation (mm)",
       x = "Week",
       fill = "Year") +
  facet_wrap(~year, scales = "free_x") +
  theme_pubr(legend = "bottom", margin = TRUE) +
  labs_pubr()

ppt_all_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_ppt.tiff",
       plot = ppt_all_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# SRER DESGR Summary Stats
describeBy(x = ppt_summary, group = c("year", "month"))

sum_2017_monsoon <- ppt_summary %>% 
  group_by(day, year) %>% 
  filter(as.integer(day) >= 166 & as.integer(day) <= 273) %>% 
  summarize(ppt_mm_sum = sum(ppt_mm),
            ppt_in_sum = sum(ppt_in)) %>% 
  filter(year == 2017) %>% 
  group_by(year) %>% 
  summarize(tot_mm = sum(ppt_mm_sum),
            tot_in = sum(ppt_in_sum))

sum_2018_monsoon <- ppt_summary %>% 
  group_by(day, year) %>% 
  filter(as.integer(day) >= 166 & as.integer(day) <= 258) %>% 
  summarize(ppt_mm_sum = sum(ppt_mm),
            ppt_in_sum = sum(ppt_in)) %>% 
  filter(year == 2018) %>% 
  group_by(year) %>% 
  summarize(tot_mm = sum(ppt_mm_sum),
            tot_in = sum(ppt_in_sum))

sum_2019_monsoon <- ppt_summary %>% 
  group_by(day, year) %>% 
  filter(as.integer(day) >= 166 & as.integer(day) <= 258) %>% 
  summarize(ppt_mm_sum = sum(ppt_mm),
            ppt_in_sum = sum(ppt_in)) %>% 
  filter(year == 2019) %>% 
  group_by(year) %>% 
  summarize(tot_mm = sum(ppt_mm_sum),
            tot_in = sum(ppt_in_sum))

monsoon_tot <- rbind(sum_2017_monsoon, sum_2018_monsoon, sum_2019_monsoon)

# Sig difference between precipitation b/w years?
hist(sqrt((log10(ppt_summary$ppt_mm))))
lm(sqrt(log10(ppt_mm))~year, data = ppt_summary)
summary(lm(sqrt(log10(ppt_mm))~year, data = ppt_summary)) # No


### SRER DESGR Air Temp
temp <- vroom("Data/site-env-data/ppt-temp/all_air_ppt_temp.csv",
                  col_types = c(.default = "d",
                                datetime = "T",
                                pptmm = "n",
                                tempC = "n",
                                ids = "f"))
str(temp)

temp_series <- temp %>% 
  dplyr::select(-pptmm) %>% 
  filter(datetime >= ymd_hms('2017-01-01 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  tidyr::extract("ids", into = c("site"), "(.{4})", remove = TRUE) %>% 
  mutate(site = as.factor(site))

glimpse(temp_series)
summary(temp_series)

temp_parse <- temp_series %>%
  drop_na(tempC) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  mutate(day = as.factor(yday(datetime)),
         week = week(datetime),
         month = as.factor(month(datetime)),
         year = as.factor(year(datetime))) %>%
  filter(year != 2020) %>% 
  group_by(site, day, week, month, year) %>% 
  arrange(site, year, month, week, day)

temp_summary <- temp_parse %>% 
  group_by(day, week, month, year) %>% 
  summarise(site_temp = mean(tempC)) %>% 
  arrange(year, month, week, day)

temp_all_fig <- temp_summary %>% 
  group_by(week, year) %>% 
  ggplot(aes(x = week, y = site_temp, color = year)) +
  geom_smooth(span = 0.07, se = FALSE) +
  #geom_line(size = 1) +
  scale_color_manual(values = c("#DB4325", "#DB4325", "#DB4325")) +
  scale_x_continuous(breaks=seq(0, 52, 4)) +
  labs(y = "Temperature (°C)",
       x = "Week",
       color = "Year") +
  facet_wrap(~year, scales = "free_x") +
  theme_pubr(legend = "none", margin = TRUE) +
  labs_pubr()

temp_all_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_temp.tiff",
       plot = temp_all_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Combo SRER DESGR PPT and Temp

site_env_all_fig <- ppt_final %>% 
  group_by(week, year) %>% 
  ggplot(aes(x = week)) +
  geom_col(aes(y = ppt_mm, fill = year)) +
  #geom_hline(yintercept = 34.25, color = "darkgreen", size = 1.5) +
  geom_smooth(mapping = aes(y = site_temp, color = "#DB4325"),
              size = 1.5, 
              span = 0.07,
              se = FALSE,
              data = temp_summary) +
  scale_color_discrete(guide = guide_legend(label = FALSE)) +
  scale_fill_manual(values = c("#b0e8f5", "#169cf0", "#0033FF")) +
  scale_x_continuous(breaks = seq(1, 52, 3),
                     expand = c(0,0)) +
  scale_y_continuous(name = "Precipitation (mm)",
                     breaks = seq(0, 130, 10),
                     sec.axis = sec_axis(~.,
                                         name = "Mean Temperature (°C)",
                                         breaks = seq(0, 40, 10)),
                     expand = c(0.01,0)) +
  labs(x = "Week",
       fill = "Precipitation",
       color = "Temperature") +
  facet_wrap(~year) +
  theme_pubr(legend = "bottom", margin = TRUE) +
  theme(axis.title.y = element_text(hjust = 0.4,
                                    vjust = 1,
                                    size = 14),
        axis.title.y.right = element_text(hjust = 0.95,
                                          vjust = 1,
                                          size = 14)) +
  labs_pubr()
  
site_env_all_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_ppt_temp.tiff",
       plot = site_env_all_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Monsoon only SRER DESGR Precipitation and Temp
site_env_monsoon_ppt <- ppt_final %>% 
  mutate(day = as.integer(day)) %>% 
  filter(day >= 166 & day <= 273)

site_env_monsoon_ppt %>% 
  group_by(year) %>% 
  summarise(sum_ppt = sum(ppt_mm))

# number of rain days >= 0.5
count_raindays <- site_env_monsoon_ppt %>% 
  group_by(day, year) %>% 
  summarise(count_ppt = sum(ppt_mm > 0))

count_raindays %>% 
  group_by(year) %>% 
  summarise(sum_days = sum(count_ppt))

count_raindays_arranged <- count_raindays %>% 
  arrange(year, day)
         
site_env_monsoon_temp <- temp_summary %>% 
  mutate(day = as.integer(day)) %>% 
  filter(day >= 166 & day <= 273)

site_env_monsoon_fig <- site_env_monsoon_ppt %>% 
  group_by(day, year) %>% 
  mutate(label = paste(month, monthday, sep = "-"),
    monthday = as.integer(monthday),
    label = as.factor(label)) %>% 
  ggplot(aes(x = label)) +
  geom_col(aes(y = ppt_mm, fill = year), color = "black") +
  # geom_hline(yintercept = 34.25, color = "darkgreen", size = 1.5) +
  # geom_smooth(mapping = aes(y = site_temp, color = "#FF3300"),
  #             size = 1.5,
  #             span = 0.07,
  #             se = FALSE,
  #             data = site_env_monsoon_temp) +
  # scale_color_discrete(guide = guide_legend(label = FALSE)) +
  # geom_smooth(mapping = aes(y = ppt_mm), color = "#330066", show.legend = TRUE,
  #             size = 1.5,
  #             span = 0.5,
  #             se = FALSE) +
  scale_fill_manual(values = c("#b0e8f5", "#0033FF", "#169cf0")) +
  scale_x_discrete(breaks = c("Jun-15", "Jun-25", "Jul-5", "Jul-15", "Jul-25",
                              "Aug-5", "Aug-15", "Aug-25", "Sep-5", "Sep-15",
                              "Sep-25")) +
  scale_y_continuous(name = "Precipitation (mm)",
                     breaks = seq(0, 150, 10),
                     expand = c(0.01,0)) +
                     # sec.axis = sec_axis(~.,
                     #                     name = "Mean Temperature (°C)",
                     #                     breaks = seq(0, 40, 10)),
                     # expand = c(0.0,0)) +
  labs(x = "Month - Day",
       fill = "Precipitation",
       color = "Temperature") +
  facet_wrap(~year, strip.position = c("bottom")) +
  theme_pubr(legend = "top", margin = T, x.text.angle = 45) +
  theme(axis.title.y = element_text(hjust = 0.4, vjust = 1)) +
  labs_pubr(base_size = 24)
        # axis.title.y.right = element_text(hjust = 0.9,
        #                                   vjust = 1,
        #                                   size = 10))


site_env_monsoon_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_monsoon_ppt.tiff",
       plot = site_env_monsoon_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

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
hist(stemp_series$dailyavg) # looks good

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
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom") +
  labs_pubr()

below_soil_temp_fig

ggsave(filename = "Figures_Tables/environment/5cm_soil_temp.tiff",
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
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom") +
  labs_pubr()

surface_temp_fig

ggsave(filename = "Figures_Tables/environment/surface_soil_temp.tiff",
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
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom") +
  labs_pubr()

light_fig

ggsave(filename = "Figures_Tables/environment/surface_light.tiff",
       plot = light_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Soil Moisture 

### NEED TO FINISH

sm <- vroom("Data/site-env-data/soil-moisture/all_sm.csv",
                    col_types = c(datetime = "T",
                                  precip = "f",
                                  clip = "f",
                                  moisture = "n",
                                  ids = "f"))
str(sm)


sm_clean <- sm %>% 
  dplyr::select(-ids) %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-18 00:00:00')) %>% 
  distinct(datetime, precip, clip, .keep_all = TRUE) %>% 
  mutate(moisture = replace(moisture, moisture < 0, NA)) %>%
  mutate(moisture = replace(moisture, moisture > 0.3, NA)) %>% 
  drop_na()

glimpse(sm_clean)
summary(sm_clean)

sm_parse <- sm_clean %>%
  mutate(day = yday(datetime),
         week = week(datetime),
         month = as.factor(month(datetime, label = T, abbr = T)),
         year = as.factor(year(datetime)),
         monthday = as.factor(mday(datetime))) %>% 
  group_by(precip, clip, month, monthday, year) %>% 
  summarise(sm_mean = mean(moisture))

sm_summary <- sm_parse %>% 
  filter(year != 2017 & year != 2020) %>% 
  group_by(precip, clip) %>% 
  summarise(sm_mean_year = mean(sm_mean))

sm_fig <- sm_parse %>% 
  group_by(precip, clip, month, monthday, year) %>% 
  filter(year != 2020) %>% 
  mutate(label = paste(month, monthday, sep = "-"),
         monthday = as.integer(monthday),
         label = as.factor(label)) %>% 
  ggplot(aes(x = label, y = 100*sm_mean, color = precip, group = precip)) +
  geom_smooth(span = 0.25, se = TRUE, size = 2) +
  geom_hline(data = sm_summary, aes(yintercept = 100*sm_mean_year, color = precip),
             size = 1, linetype = 2) +
  scale_color_manual(values = c("grey30", "blue1", "#ba7525"),
                     labels = c("Ambient", "Wet", "Drought")) +
  scale_x_discrete(breaks = c("Jan-1", "Feb-1", "Mar-1", "Apr-1", "May-1",
                              "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1",
                              "Nov-1", "Dec-1")) +
  labs(y = "Volumetric Water Content (%)",
       x = "Month - Day",
       color = "Precipitation") +
  facet_wrap(~year, scales = "free_x", nrow=1) +
  theme_pubr(legend = "right", margin = TRUE, x.text.angle = 45) +
  labs_pubr(base_size = 24)

sm_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_sm_pptx.tiff",
       plot = sm_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")



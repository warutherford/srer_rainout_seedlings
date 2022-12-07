# USDA Rainout - Santa Rita Experimental Range
# Data visualizations-site environmental data
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-05-10

# Load packages
library(tidyverse)
library(ggbreak)
library(agricolae)
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

ppt_srer_avg_all <- ppt_srer %>%
  group_by(YEAR) %>% 
  rowwise() %>% 
  dplyr::summarise(ppt_sum = sum(c_across(JAN:DEC), na.rm = TRUE)) %>% 
  mutate(ppt_yr_in = ppt_sum/100) %>% 
  mutate(ppt_yr_mm = ppt_yr_in*25.4)

summary(ppt_srer_avg_all)

ppt_srer_avg_monsoon <- ppt_srer %>%
  group_by(YEAR) %>% 
  rowwise() %>% 
  dplyr::summarise(ppt_sum = sum(c_across(JUN:SEP), na.rm = TRUE)) %>% 
  mutate(ppt_yr_in = ppt_sum/100) %>% 
  mutate(ppt_yr_mm = ppt_yr_in*25.4)

summary(ppt_srer_avg_monsoon)


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

# stats
summary(glmmTMB(ppt_mm ~ year,
        data = ppt_final,
        ziformula = ~.,
        family = poisson(link = "log")))

ppt_clean <- ppt_final %>% filter(ppt_mm != 0)
hist(log(ppt_clean$ppt_mm))
ppt.aov <- aov(log(ppt_mm)~year, data = ppt_clean)
summary(ppt.aov)
ppt.aov.post <- HSD.test(ppt.aov, "year", group = FALSE, console = TRUE)


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
  filter(as.integer(day) >= 166 & as.integer(day) <= 273) %>% 
  summarize(ppt_mm_sum = sum(ppt_mm),
            ppt_in_sum = sum(ppt_in)) %>% 
  filter(year == 2018) %>% 
  group_by(year) %>% 
  summarize(tot_mm = sum(ppt_mm_sum),
            tot_in = sum(ppt_in_sum))

sum_2019_monsoon <- ppt_summary %>% 
  group_by(day, year) %>% 
  filter(as.integer(day) >= 166 & as.integer(day) <= 273) %>% 
  summarize(ppt_mm_sum = sum(ppt_mm),
            ppt_in_sum = sum(ppt_in)) %>% 
  filter(year == 2019) %>% 
  group_by(year) %>% 
  summarize(tot_mm = sum(ppt_mm_sum),
            tot_in = sum(ppt_in_sum))

monsoon_tot <- rbind(sum_2017_monsoon, sum_2018_monsoon, sum_2019_monsoon)


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
         month = as.factor(month(datetime, label = TRUE)),
         monthday = as.factor(mday(datetime)),
         year = as.factor(year(datetime))) %>%
  filter(year != 2020) %>% 
  group_by(site, day, week, month, monthday, year) %>% 
  arrange(site, year, month, week,monthday, day)

temp_summary <- temp_parse %>% 
  group_by(day, week, monthday, month, year) %>% 
  summarise(site_temp = mean(tempC)) %>% 
  arrange(year,week,day)

temp_min_day <- temp_summary %>% 
  group_by(day, month, year) %>%
  summarise(site_monthday_temp = min(site_temp)) %>%
  arrange(year, month, day) %>% 
  group_by(day, year) %>% 
  summarise(site_temp_min = min(site_monthday_temp))

temp_year <- temp_summary %>% 
  group_by(year) %>% 
  summarise(site_temp_mean = mean(site_temp))

temp_min_yr <- temp_summary %>% 
  group_by(day, year) %>%
  summarise(site_month_temp = mean(site_temp)) %>% 
  group_by(year) %>% 
  summarise(site_temp_min = min(site_month_temp))

temp_max_yr <- temp_summary %>% 
  group_by(day, year) %>%
  summarise(site_month_temp = mean(site_temp)) %>% 
  group_by(year) %>% 
  summarise(site_temp_max = max(site_month_temp))

temp_summary %>% 
  ggplot(aes(x = day, y = site_temp, group = year, color = year)) + 
  geom_line()+
  facet_wrap(~year)+
  theme_pubr()

temp_monsoon <- temp_parse %>% 
  group_by(month, year) %>% 
  summarise(site_temp = mean(tempC)) %>% 
  arrange(year,month) %>% 
  filter(month == "Jun"| month =='Jul'|month=='Aug'|month=='Sep') %>% 
  group_by(year) %>% 
  summarise(site_temp_mon <- mean(site_temp))

temp_min_monsoon <- temp_parse %>% 
  group_by(month, year) %>% 
  summarise(site_temp = mean(tempC)) %>% 
  arrange(year,month) %>% 
  filter(month == "Jun"| month =='Jul'|month=='Aug'|month=='Sep') %>% 
  group_by(year) %>% 
  summarise(site_temp_min = min(site_temp))

temp_max_monsoon <- temp_parse %>% 
  group_by(month, year) %>% 
  summarise(site_temp = mean(tempC)) %>% 
  arrange(year,month) %>% 
  filter(month == "Jun"| month =='Jul'|month=='Aug'|month=='Sep') %>% 
  group_by(year) %>% 
  summarise(site_temp_max = max(site_temp))

# long term temp from prism (1932-Feb 2020)
temp_prism <- vroom("Data/site-env-data/prism_desgr_temp.csv",
              col_types = c(.default = "f",
                            date = NULL,
                            year = "f",
                            month = "f",
                            tmin = "n",
                            tmean = "n",
                            tmax = "n"))

temp_year_prism <- temp_prism %>% 
  group_by(year) %>% 
  summarise(site_tmin_mean = mean(tmin),
            site_tmean_mean = mean(tmean),
            site_tmax_mean = mean(tmax))

min<-temp_prism %>% group_by(month,year) %>% summarise(min = min(tmin))

long_term <- temp_year_prism %>% 
  summarise(tmin_mean = mean(site_tmin_mean),
         tmean_mean = mean(site_tmean_mean),
         tmax_mean = mean(site_tmax_mean))

temp_monsoon <- temp_prism %>% 
  group_by(month, year) %>% 
  filter(month == "6"| month =='7'|month=='8'|month=='9') %>% 
  group_by(year) %>% 
 summarise(site_tmin_mean = mean(tmin),
                      site_tmean_mean = mean(tmean),
                      site_tmax_mean = mean(tmax))

long_term_monsoon <- temp_monsoon %>% 
  summarise(tmin_mean = mean(site_tmin_mean),
            tmean_mean = mean(site_tmean_mean),
            tmax_mean = mean(site_tmax_mean))


# stats
hist(temp_summary$site_temp)
temp.aov <- aov(site_temp~year, data = temp_summary)
summary(temp.aov)
temp.aov.post <- HSD.test(temp.aov, "year", group = FALSE, console = TRUE)
# 2017 was statistically hotter than 2019, but not 2018

### Combo SRER DESGR PPT and Temp

ppt_lab <- ppt_final %>% 
  group_by(week, year) %>% 
  mutate(label = paste(month, monthday, sep = "-"),
         monthday = as.integer(monthday),
         label = as.factor(label))

temp_lab <- temp_summary %>% 
  group_by(week, year) %>% 
  mutate(label = paste(month, monthday, sep = "-"),
         monthday = as.integer(monthday),
         label = as.factor(label))

all_env <- full_join(ppt_lab, temp_lab, by = c("label", "year"))

site_env_all_fig <- all_env %>% 
  #group_by(week.x, year) %>% 
  ggplot(aes(x = week.x)) +
  geom_col(aes(y = ppt_mm, fill = year)) +
  #geom_hline(yintercept = 34.25, color = "darkgreen", size = 1.5) +
  geom_smooth(mapping = aes(y = site_temp, color = "#DB4325"),
              size = 1.5, 
              span = 0.07,
              se = FALSE) +
  scale_color_discrete(guide = guide_legend(label = FALSE)) +
  scale_fill_manual(values = c("#b0e8f5", "#169cf0", "#0033FF")) +
  scale_x_continuous(breaks = c(1, 10, 18, 27, 36, 45, 52), label = c("Jan-1", "Mar-1",
                                                                  "May-1", "Jul-1",
                                                                  "Sep-1","Nov-1", "Dec-31")) +
  scale_y_continuous(name = "Precipitation (mm)",
                     breaks = seq(0, 130, 10),
                     sec.axis = sec_axis(~.,
                                         name = "Mean Temperature (°C)",
                                         breaks = seq(0, 30, 10)),
                     expand = c(0.01,0)) +
  labs(x = "Month - Day",
       fill = "Precipitation",
       color = "Temperature") +
  facet_wrap(~year, strip.position = c("bottom")) +
  theme_pubr(legend = "bottom", margin = TRUE, x.text.angle = 45) +
  theme(axis.title.y = element_text(hjust = 0.4,
                                    vjust = 1,
                                    size = 24),
        axis.title.y.right = element_text(hjust = 1,
                                          vjust = 1,
                                          size = 24)) +
  labs_pubr(base_size = 24)
  
site_env_all_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_ppt_temp_figs1.tiff",
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

# stats, zero-inflated
summary(glmmTMB(ppt_mm ~ year,
                data = site_env_monsoon_ppt,
                ziformula = ~.,
                family = poisson(link = "log")))

hist(log(site_env_monsoon_ppt$ppt_mm))
monsoon.aov <- aov(ppt_mm~year, data = site_env_monsoon_ppt)
summary(monsoon.aov)
monsoon.aov.post <- HSD.test(monsoon.aov, "year", group = FALSE, console = TRUE)

# summary rain calcs
site_env_monsoon_ppt %>% 
  group_by(year) %>% 
  summarise(sum_ppt = sum(ppt_mm))

# number of rain days >= 0.5
count_raindays <- site_env_monsoon_ppt %>% 
  group_by(day, year) %>% 
  summarise(count_ppt = sum(ppt_mm >= 0.5))

count_raindays %>% 
  group_by(year) %>% 
  summarise(sum_days = sum(count_ppt))

count_raindays_arranged <- count_raindays %>% 
  arrange(year, day)

# get only 2017 values
consec_raindays_2017 <- count_raindays_arranged %>% ungroup() %>% 
  filter(year == 2017) %>% dplyr::select(-year,-day)

# highest number with value 1, which indicates a day with rain >= 0.5 mm
consec_seq_2017 <- data.frame(unclass(rle(consec_raindays_2017$count_ppt)))

# largest gap = 9 days not counting the end of the monsoon dry down

# get only 2018 values
consec_raindays_2018 <- count_raindays_arranged %>% ungroup() %>% 
  filter(year == 2018) %>% dplyr::select(-year,-day)

# highest number with value 1, which indicates a day with rain >= 0.5 mm
consec_seq_2018 <- data.frame(unclass(rle(consec_raindays_2018$count_ppt)))  

# largest gap = 21 days between rain events >= 0.5mm, followed by 14

# get only 2019 values
consec_raindays_2019 <- count_raindays_arranged %>% ungroup() %>% 
  filter(year == 2019) %>% dplyr::select(-year,-day)

# highest number with value 1, which indicates a day with rain >= 0.5 mm
consec_seq_2019 <- data.frame(unclass(rle(consec_raindays_2019$count_ppt))) 

# largest gap = 11 days not counting the dry start of the monsoon

# 2017 = 16; 2018 = 8; 2019 = 4 consecutive days with rain (>= 0.5 mm)

## site monsoon ppt for each year
site_env_monsoon_temp <- temp_summary %>% 
  mutate(day = as.integer(day)) %>% 
  filter(day >= 166 & day <= 273)

lines_dat <- data.frame(year = c(2017,2018, 2019), precip = c(280, 330, 292)) %>% 
  mutate(year = as.factor(year))

site_env_monsoon_fig <- site_env_monsoon_ppt %>% 
  group_by(day, year) %>% 
  mutate(label = paste(month, monthday, sep = "-"),
    monthday = as.integer(monthday),
    label = as.factor(label)) %>% 
  ggplot(aes(x = label)) +
  geom_col(aes(y = ppt_mm, fill = year), color = "black", show.legend = FALSE) +
  geom_hline(data = lines_dat, aes(yintercept = precip, group = year, color = year), size = 2.5) +
  geom_hline(aes(yintercept = 240, color = "Long-term Mean"), size = 2.5, show.legend = TRUE) +
  scale_fill_manual(values = c("#b0e8f5", "#0033FF", "#169cf0"), na.value = "") +
  scale_color_manual(values = c("#b0e8f5", "#0033FF", "#169cf0", "darkorange")) +
  scale_x_discrete(breaks = c("Jun-15", "Jun-30", "Jul-15", "Jul-30",
                             "Aug-15", "Aug-31", "Sep-15", "Sep-30")) +
  ylim(0, 350)+
  scale_y_break(breaks = c(70, 100),scales = "free", ticklabels = c(100, 105, 110), space = 0.5)+
  scale_y_break(breaks = c(110, 200), scales = "free", space = 0.5)+
  labs(x = "Month - Day",
       y = "PPT (mm)",
       fill = NULL,
       color = "Total Monsoon PPT",
       yintercept = NULL) +
  facet_wrap(~year, strip.position = c("bottom")) +
  theme_pubr(legend = "top", margin = T, x.text.angle = 45) +
  theme(axis.title.y = element_text(vjust = -1), panel.spacing.x = unit(0.5, "inches")) +
  labs_pubr(base_size = 24)+
  guides(fill = guide_legend(override.aes = list(linetype = 0), label = FALSE))

site_env_monsoon_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_monsoon_ppt_fig1.tiff",
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
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  group_by(precip, location, clip) %>% 
  summarise(tempC_mean = mean(soiltempC),
            tempC_se = sd(soiltempC)/sqrt(n()))

# make sure date range is within seedling observation date range
stemp_stats <- stemplight %>% 
  dplyr::select(-lux) %>% 
  filter(location == 'BelowTemp') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE)

# pull out only surface soil temp
atemp <- stemplight %>% 
  dplyr::select(-lux) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00')) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  group_by(precip, location, clip) %>% 
  summarise(atempC_mean = mean(soiltempC),
            atempC_se = sd(soiltempC)/sqrt(n()))

# time series of 5cm soil temp
stemp_series <- stemplight %>% 
  dplyr::select(-lux) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  filter(location == 'BelowTemp') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00')) %>% 
  mutate(hour = hour(datetime),
         week = isoweek(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>%
  summarise(dailyavg = mean(soiltempC))

# test sig diff of daily below ground soil temp between treatments
# data normal?
hist(stemp_series$dailyavg) # looks good

stemp_series %>% group_by(precip, year) %>% summarise(meant = mean(dailyavg))

below_anova <- aov(dailyavg ~ precip*clip*year, data = stemp_series)

summary(below_anova)

below.aov.post <- HSD.test(below_anova, "year", group = FALSE, console = TRUE)

# figure of 5-cm soil temp time series
below_soil_temp_fig <- stemp_series %>% 
  group_by(precip, clip, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg, color = clip))+
  geom_line(size = 2)+
  scale_x_continuous(breaks = c(1, 10, 18, 27, 36, 45, 52), label = c("Jan-1", "Mar-1",
                                                                      "May-1", "Jul-1",
                                                                      "Sep-1","Nov-1", "Dec-31")) +
  scale_color_manual(values = c("brown","darkblue")) +
  #xlim(0,52)+
  labs(y = "Daily 5-cm Soil Temperature (°C)",
       x = "Month-Day",
       color = "Grazing/Clipping",
       linetype = "Year") +
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom", margin = F, x.text.angle = 45) +
  labs_pubr(base_size = 24)

below_soil_temp_fig

ggsave(filename = "Figures_Tables/environment/5cm_soil_temp_figs4.tiff",
       plot = below_soil_temp_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Soil Surface Temp

# time series of surface soil temp, make sure time is within seedling observations
atemp_series <- stemplight %>% 
  dplyr::select(-lux) %>% 
  distinct(datetime, .keep_all = TRUE) %>% 
  filter(location == 'TempLight') %>% 
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00'))

# Corrupted sensor data make NA
atemp_fix <- atemp_series %>% 
  filter(precip == 'Ambient' & section == 'E1') %>% 
  filter((datetime >= ymd('2019-05-06') & datetime <= ymd('2019-08-01'))) %>% 
  mutate(soiltempC = replace(soiltempC, soiltempC > 0, NA))

# make a dummy dataset of missing sensor data to keep lines from connecting on soil surface temp and light figs
gap_fix <- vroom("Data/site-env-data/soil-temp/gap_fix.csv",
                 col_types = c(.default = "d",
                               precip = "f",
                               clip = "f",
                               week = "d",
                               year = "f",
                               dailyavg = "d"))

# join full data set with NAs, extract daily averages for graphing
atemp_comp <- full_join(atemp_series, atemp_fix) %>% 
  mutate(hour = hour(datetime),
         week = isoweek(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>% 
  summarise(dailyavg = mean(soiltempC))

atemp_comp <- full_join(atemp_comp, gap_fix)

# test sig diff of daily surface soil temp between treatments
# data normal?
hist(atemp_comp$dailyavg) # normal
atemp_comp %>% group_by(precip, year) %>% summarise(meant = mean(dailyavg, na.rm = T))


atemp_aov <- aov(dailyavg ~ precip*clip*year, data = atemp_comp)

summary(atemp_aov)

post.hoc.surf <- HSD.test(atemp_aov, "precip", group = FALSE, console = TRUE)
post.hoc.surf

# plot time series with missing data in ambient clipped from 2019
surface_temp_fig <- atemp_comp %>%
  group_by(precip, clip, week, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg, color = clip)) +
  geom_line(size = 2) +
  scale_x_continuous(breaks = c(1, 10, 18, 27, 36, 45, 52), label = c("Jan-1", "Mar-1",
                                                                      "May-1", "Jul-1",
                                                                      "Sep-1","Nov-1", "Dec-31")) +
  scale_color_manual(values = c("brown","darkblue")) +
  #xlim(NA, 52) +
  labs(y = "Daily Soil Suface Temperature (°C)",
       x = "Month - Day",
       color = "Grazing/Clipping",
       linetype = "Year") +
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom", margin = F, x.text.angle = 45) +
  labs_pubr(base_size = 24)

surface_temp_fig

ggsave(filename = "Figures_Tables/environment/surface_soil_temp_figs3.tiff",
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
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00') & 
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
  filter(datetime >= ymd_hms('2017-07-15 00:00:00') & datetime <= ymd_hms('2020-01-01 00:00:00') & 
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
         week = isoweek(datetime),
         year = as.factor(year(datetime))) %>% 
  group_by(precip, clip, week, year) %>% 
  summarise(dailyavg_light = mean(lux))

light_comp <- full_join(light_comp, gap_fix)

# test sig diff of daily below ground soil temp between treatments
# data normal?
hist(light_comp$dailyavg_light) # looks good

light_anova <- aov(dailyavg_light ~ precip*clip*year, data = light_comp)

summary(light_anova)

post.hoc.light <- HSD.test(light_anova, "precip", group = FALSE, console = TRUE)
post.hoc.light

# figure of light time series
light_fig <- light_comp %>% 
  group_by(precip, clip, year) %>% 
  ggplot(mapping = aes(x = week, y = dailyavg_light, color = clip))+
  geom_line(size = 2)+
  scale_x_continuous(breaks = c(1, 10, 18, 27, 36, 45, 52), label = c("Jan-1", "Mar-1",
                                                                      "May-1", "Jul-1",
                                                                      "Sep-1","Nov-1", "Dec-31")) +
  scale_y_continuous(limits = c(0, 180000), labels = scales::scientific) + 
  scale_color_manual(values = c("brown","darkblue")) +
  #ylim(0,150000)+
  labs(y = "Daily Light (lux)",
       x = "Month - Day",
       color = "Grazing/Clipping",
       linetype = "Year") +
  facet_wrap(~precip+year, scales = "free_x") +
  theme_pubr(legend = "bottom", margin = F, x.text.angle = 45) +
  labs_pubr(base_size = 24)

light_fig

ggsave(filename = "Figures_Tables/environment/surface_light_figs5.tiff",
       plot = light_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Soil Moisture 

sm_1 <- vroom("Data/site-env-data/soil-moisture/all_sm_1.csv",
                    col_types = c(datetime = "T",
                                  precip = "f",
                                  clip = "f",
                                  moisture = "n",
                                  ids = "f"))

sm_2 <- vroom("Data/site-env-data/soil-moisture/all_sm_2.csv",
              col_types = c(datetime = "T",
                            precip = "f",
                            clip = "f",
                            moisture = "n",
                            ids = "f"))

sm <- bind_rows(sm_1, sm_2)

str(sm)

sm_clean <- sm %>% 
  dplyr::select(-ids) %>% 
  filter(datetime >= ymd_hms('2018-01-01 00:00:00') & datetime <= ymd_hms('2019-12-31 00:00:00')) %>% 
  distinct(datetime, precip, clip, .keep_all = T) %>% 
  mutate(moisture = replace(moisture, moisture < 0, NA)) %>%
  mutate(moisture = replace(moisture, moisture > 0.3, NA)) %>% 
  drop_na()

glimpse(sm_clean)
summary(sm_clean)

# stats
hist(sm_clean$moisture)
sm.aov <- aov(moisture~precip*clip, data = sm_clean)
summary(sm.aov)
sm.aov.post <- HSD.test(sm.aov, "precip", group = F, console = TRUE)

# make sm fig
sm_parse <- sm_clean %>%
  mutate(day = yday(datetime),
         week = week(datetime),
         month = as.factor(month(datetime, label = T, abbr = T)),
         year = as.factor(year(datetime)),
         monthday = as.factor(mday(datetime))) %>% 
  group_by(precip, month, monthday, year) %>% 
  summarise(sm_mean = mean(moisture))

sm_good <- sm_parse %>%  filter(precip != "IR" | year != "2018")

sm_replace <- sm_parse %>% 
  filter(precip == "IR" & year == "2018" & month > "May")

sm_fix <- rbind(sm_good, sm_replace)

sm_summary <- sm_fix %>% 
  group_by(precip, year) %>% 
  summarise(sm_mean_year = mean(sm_mean))

sm_fig <- sm_fix %>% 
  group_by(precip, month, monthday, year) %>% 
  mutate(label = paste(month, monthday, sep = "-"),
         monthday = as.integer(monthday),
         label = as.factor(label)) %>% 
  ggplot(aes(x = label, y = 100*sm_mean, color = precip, group = precip)) +
  geom_smooth(span = 0.3, se = F, size = 2) +
  geom_hline(data = sm_summary, aes(yintercept = 100*sm_mean_year, color = precip),
             size = 2, linetype = 2) +
  scale_color_manual(values = c("grey30", "blue1", "#ba7525"),
                     labels = c("Ambient", "Wet", "Drought")) +
  scale_x_discrete(breaks = c("Jan-1", "Feb-1", "Mar-1", "Apr-1", "May-1",
                              "Jun-1", "Jul-1", "Aug-1", "Sep-1", "Oct-1",
                              "Nov-1", "Dec-1")) +
  labs(y = "Volumetric Water Content (%)",
       x = "Month - Day",
       color = "PPTx") +
  facet_wrap(~year, ncol = 2, strip.position = "bottom") +
  theme_pubr(legend = "right", margin = TRUE, x.text.angle = 45) +
  labs_pubr(base_size = 24)

sm_fig

ggsave(filename = "Figures_Tables/environment/srer_desgr_sm_pptx_figs2.tiff",
       plot = sm_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

### Rainout PAR measurements with line quantum sensor
par <- vroom("Data/site-env-data/par/rainout-lqs-complete.csv",
                  col_types = c(.default = "d",
                                date = "D",
                                block = "f",
                                tx = "f",
                                location = "f",
                                clip = "f"))
str(par)

# convert from wide to long
par_long <- par %>% 
  pivot_longer(6:8, names_to = "transect", values_to = "value")

# calculate average of three transects
par_calc <- par_long %>% 
  group_by(date, block, tx, location, clip) %>% 
  summarise(mean_read = mean(value))

# separate above reads with below canopy ready to create ratio of intercepted light
par_calc_above <- par_calc %>% 
  filter(location == "Above") %>% 
  droplevels()

par_calc_below <- par_calc %>% 
  filter(location == "Below" & tx != "RO-Above") %>% 
  droplevels()

par_loc <- bind_cols(par_calc_above, par_calc_below, .name_repair = "unique") %>% 
  dplyr::select(1:6,10,12) %>% 
  rename(date = date...1,
         block = block...2,
         precip = tx...3,
         above = location...4,
         clip = clip...5,
         abv_read = mean_read...6,
         below = location...10,
         bel_read = mean_read...12)

# calculate intercepted PAR ratio
par_ratio <- par_loc %>%  mutate(intercept_ratio = 1-(bel_read/abv_read))

summary(par_ratio)

# look at summary of ratios
par_ratio %>% group_by(precip, clip) %>% 
  summarise(mean_inter = mean(intercept_ratio))

# create violin plot of distribution of measurements by ppt and clip
par_fig <- par_ratio %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = clip, y = 100*intercept_ratio, fill = clip)) +
  geom_violin() +
  #geom_boxplot()+
  scale_fill_manual(values = c("darkorange","forestgreen")) +
  scale_x_discrete(labels = c("","")) +
  labs(y = "Intercepted PAR (%)",
       x = "") +
  ylim(NA,100)+
  theme_pubr(legend = "none") +
  facet_wrap(~precip)+
  labs_pubr(base_size = 24)

par_fig

# save plot
ggsave(filename = "Figures_Tables/environment/lqs_par_figs6.tiff",
       plot = par_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# Anova of differences between precip and clipping txs
hist(sqrt(par_ratio$intercept_ratio)) # better with sqrt

par_aov <- aov(sqrt(intercept_ratio)~precip*clip, data = par_ratio)

summary(par_aov)

post.par <- HSD.test(par_aov, c("precip", "clip"), group = T, console = TRUE)



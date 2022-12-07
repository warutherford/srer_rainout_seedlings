# USDA Rainout - Santa Rita Experimental Range
# Stats and Modeling
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-09-29

# Load packages
library(tidyverse)
library(vroom)
library(psych)
library(mgcv)
library(lmerTest)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(multcomp)

# The data cleaning and formatting code below is also at the top of each response
# model script as well. This script covers how the final data set was prepped to do
# modeling, followed by the descriptive stats
# Herbivory not accessed as part of this study.

### Read in seedlings data (cohorts 1-3), make all columns factor and date a date
seedlings <- vroom("Data/seedlings_combined.csv",
                   col_select = -c(1),
                   col_types = c(.default = "f",
                                 date = "D"))
str(seedlings)

# Fate counts
seedlings_fate_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(fate) %>% 
  pivot_wider(names_from = fate,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_germ = "0",
         survival = "1",
         died = "2") %>% 
  mutate(tot_germination = survival + died,
         surv_perc = (survival/tot_germination)) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0))

hist((seedlings_fate_full$surv_perc))

seedling_fate_post_co1y1 <- seedlings_fate_full %>% 
  filter(cohort == 1) %>% 
  filter(date > "2017-07-21" & date < "2018-07-10") %>%  # emerge (2 weeks following planting) to planting of cohort 2
  mutate(year = as.factor(1))

cohort1_year1 <- seedling_fate_post_co1y1 %>% 
  group_by(cohort, precip, clip, excl, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv) # get percentage survival mean out of number that germ 
  
seedling_fate_post_co1y2 <- seedlings_fate_full %>% 
  filter(cohort == 1) %>% 
  filter(date > "2018-07-10" & date < "2019-08-01") %>%  # planting of cohort 2 to planting of cohort 3
  mutate(year = as.factor(2))

cohort1_year2 <-seedling_fate_post_co1y2 %>% 
  group_by(cohort, precip, clip, excl, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)  # get percentage survival mean out of number that germ 
  
seedling_fate_post_co1y3 <- seedlings_fate_full %>% 
  filter(cohort == 1) %>% 
  filter(date > "2019-08-01") %>% # planting of cohort 3 to end
  mutate(year = as.factor(3))

cohort1_year3 <-seedling_fate_post_co1y3 %>% 
  group_by(cohort, precip, clip, excl, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)   # get percentage survival mean out of number that germ 

seedling_fate_post_co2y1 <- seedlings_fate_full %>% 
  filter(cohort == 2) %>% 
  filter(date > "2018-07-26" & date < "2019-08-01") %>%  # emerge to planting of cohort 3
  mutate(year = as.factor(1))

cohort2_year1 <-seedling_fate_post_co2y1 %>% 
  group_by(cohort, precip, clip, excl, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv) # get percentage survival mean out of number that germ 
  
seedling_fate_post_co2y2 <- seedlings_fate_full %>% 
  filter(cohort == 2) %>% 
  filter(date > "2019-08-01") %>%  # planting of cohort three to end
  mutate(year = as.factor(2))

cohort2_year2 <-seedling_fate_post_co2y2 %>% 
  group_by(cohort, precip, clip, excl, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv) # get percentage survival mean out of number that germ 
  
seedling_fate_post_co3y1 <- seedlings_fate_full %>% 
  filter(cohort == 3) %>% 
  filter(date > "2019-08-17") %>%  # emerge of cohort 3 to end
  mutate(year = as.factor(1))

cohort3_year1 <-seedling_fate_post_co3y1 %>% 
  group_by(cohort, precip, clip, excl, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv) # get percentage survival mean out of number that germ 

surv_cohort_year <- rbind(cohort1_year1, cohort1_year2, cohort1_year3, cohort2_year1,
                          cohort2_year2, cohort3_year1)

seedling_fate_year <- rbind(seedling_fate_post_co1y1, seedling_fate_post_co1y2,
                            seedling_fate_post_co1y3, seedling_fate_post_co2y1,
                            seedling_fate_post_co2y2, seedling_fate_post_co3y1)

# Granivory counts
seedlings_gran_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(granivory) %>% 
  pivot_wider(names_from = granivory,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_gran = "0",
         granivory = "1")

# Joining fate, herbivory, and granivory variables of interest to one table
seedlings_all_full <- seedlings_fate_full %>% 
  left_join(seedlings_herb_full) %>% 
  left_join(seedlings_gran_full) %>% 
  dplyr::select(-c(no_germ, no_herb, no_gran))

# Histogram of variables, all zero-inflated poisson except for tot_germination
seedlings_all_full %>% 
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins = 50, binwidth = 1)

glimpse(seedlings_all_full)
summary(seedlings_all_full)

# descriptive stats for each cohort 1-3
describeBy(seedlings_all_full, group = "cohort")

# make a new ID column for each observation and plot to test random factors
seedlings_obs <- seedlings_all_full %>% 
  rowid_to_column("ObsID") %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:excl, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         ObsID = as.factor(ObsID),
         sampID = as.factor(sampID))

seedlings_obs_year <- seedling_fate_year %>% 
  rowid_to_column("ObsID") %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:excl, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         ObsID = as.factor(ObsID),
         sampID = as.factor(sampID))

#write.csv(seedlings_obs_year, file = "Data/seedlings_obs_year.csv")

# mixed effects model with nesting, date (cohort) within year for temporal autocorrelation and sampling as random

hist(log(seedlings_obs$survival)) # still zero-inflated
hist(sqrt(seedlings_obs$survival)) # still zero-inflated

# count data...poisson and zero-inflated (41% of data are zeros)
num_obs <- seedlings_obs %>% ungroup() %>% summarise(obs = n())

zeros <- seedlings_obs %>%
  ungroup() %>% 
  filter(survival == "0") %>% 
  summarise(zero = 100*(n()/num_obs$obs))

# glmmTMB function to build zero-inflated model, poisson dist., random factors
# see specific response variable stat scripts for specific models

# pull out descriptive stats for each treatment individually and then all treatments
# means and standard errors, not in percentages
describeBy(seedlings_obs~precip, mat = T, digits = 4)

describeBy(seedlings_obs~clip, mat = T, digits = 4)

describeBy(seedlings_obs~excl, mat = T, digits = 4)

describeBy(seedlings_obs_germ~precip+clip+excl+cohort, mat = T, digits = 4) #seedlings_obs_germ from "2_srer_rainout_model_stats_germ.R"

describeBy(seedling_fate_year~precip+clip+excl+cohort, mat = T, digits = 4)

describeBy(seedling_fate_year~clip + year, mat = T, digits = 4)

describeBy(seedling_fate_year~precip+clip+excl+cohort+year, mat = T, digits = 4)

## get confidence intervals for each treatment and then all treatments
# germination in percent
Rmisc::group.CI(100*(tot_germination/10) ~ precip,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(tot_germination/10) ~ clip,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(tot_germination/10) ~ excl,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(tot_germination/10) ~ cohort,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(tot_germination) ~ precip+clip+excl+cohort,
                data = seedlings_obs_germ,
                ci = 0.95)

# survival in percent
Rmisc::group.CI(100*(survival/10) ~ precip+clip+excl,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(survival/10) ~ precip,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(survival/10) ~ clip,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(survival/10) ~ cohort,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(survival/10) ~ excl,
                data = seedlings_obs,
                ci = 0.95)

# post hoc survival, change filters for exact tx combo for table s3
clean <- seedling_fate_year %>% filter(year != "3") %>% filter(year == "1" & precip == "RO"&clip=="Unclipped")
sur_mod <- aov(surv_perc~excl*cohort, data = clean)
TukeyHSD(sur_mod)



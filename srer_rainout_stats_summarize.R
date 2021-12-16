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
  mutate(tot_germination = survival + died)

# Herbivory counts
seedlings_herb_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(herbivory) %>% 
  pivot_wider(names_from = herbivory,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_herb = "0",
         herb_lived = "1",
         herb_died = "2") %>% 
  mutate(tot_herbivory = herb_lived + herb_died)

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
  keep(is.integer) %>% 
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

describeBy(seedlings_obs~precip+clip+excl, mat = T, digits = 4)


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

Rmisc::group.CI(100*(tot_germination/10) ~ precip+clip+excl,
                data = seedlings_obs,
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

Rmisc::group.CI(100*(survival/10) ~ excl,
                data = seedlings_obs,
                ci = 0.95)


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

# figure of mean survival and cohort by year of survival

surv_year_fig <- surv_cohort_year  %>% 
  group_by(cohort, precip, year) %>% 
  summarize(survival = mean(mean_surv)) %>% 
  mutate(precip = recode_factor(precip,
                              "Control" = "Ambient",
                              "RO" = "Drought",
                              "IR" = "Wet")) %>%
  mutate(cohort = recode_factor(cohort,
                                "1" = "2017",
                                "2" = "2018",
                                "3" = "2019")) %>%
  # mutate(excl = recode_factor(excl, 
  #                             "Control" = "None",
  #                             "Ants" = "Ants Excl",
  #                             "Rodents" = "Small Mammals Excl",
  #                             "Total" = "All Excl")) %>% 
  ggplot(aes(x = as.integer(year), y = 100*survival, color = cohort)) + 
  # geom_pointrange(aes(ymin = 10*lower, ymax = 10*upper, color = excl), size = 0.5) +
  geom_smooth(method = "glm", formula = y ~ log(x) + x, se = F, size = 2) +
  geom_point(size = 10) +
  scale_color_manual(values = c("#b0e8f5", "#169cf0", "#0033FF")) +
  labs(y = "Seedling Survival (%)",
       x = "Year of Survival",
       color = "Cohort") +
  #ylim(0, 50) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_pubr(legend = "right") +
  facet_wrap(~precip) +
  theme(panel.spacing.x = unit(2, "lines")) +
  labs_pubr(base_size = 30)

surv_year_fig

ggsave(filename = "Figures_Tables/seedlings/survival_year.tiff",
       plot = surv_year_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

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

describeBy(seedlings_obs~cohort + precip, mat = T, digits = 4)

t<-describeBy(seedling_fate_year~precip+clip+excl+cohort+year, mat = T, digits = 4)

describeBy(seedling_fate_year~precip+cohort+year, mat = T, digits = 4)

describeBy(seedling_fate_year~precip+clip+cohort+year, mat = T, digits = 4)

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

Rmisc::group.CI(100*(survival/10) ~ cohort,
                data = seedlings_obs,
                ci = 0.95)

Rmisc::group.CI(100*(survival/10) ~ excl,
                data = seedlings_obs,
                ci = 0.95)


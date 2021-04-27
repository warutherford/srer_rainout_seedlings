# USDA Rainout - Santa Rita Experimental Range
# Data visualizations
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-03-29

# Load packages
library(tidyverse)
library(forcats)
library(vroom)
library(ggpubr)
library(glmmTMB)

### Read in clean seedlings data with model outputs from "srer_rainout_stats_summarize"
seedlings_obs <- vroom("Data/seedlings_obs.csv",
                   col_types = c(.default = "f",
                                 ObsID = "i",
                                 date = "D",
                                 survival = "i",
                                 died = "i",
                                 tot_germination = "i",
                                 herb_lived = "i",
                                 herb_died = "i",
                                 tot_herbivory = "i",
                                 granivory = "i",
                                 res_surv = "d",
                                 sim_surv = "d",
                                 sim_fit_surv = "d",
                                 pred_surv = "d"))
str(seedlings_obs)
glimpse(seedlings_obs)

# create data set for raw survival data
tot_surv <- seedlings_obs %>% 
  group_by(precip, cohort, date) %>% 
  summarise(mean_surv = mean(survival),
            se_surv = (sd(survival)/n())) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv,
         date = as.Date(date))

summary(tot_surv)

tot_surv %>%
  group_by(precip, cohort) %>% 
  summarise(mean_precip_per = (100*mean(mean_surv)/10),
            mean_se_per = (100*mean(se_surv)/10)) %>% 
  mutate(upper = mean_precip_per + mean_se_per,
         lower = mean_precip_per - mean_se_per)

tot_surv %>%
  group_by(precip) %>% 
  summarise(mean_precip_per = (100*mean(mean_surv)/10),
            mean_se_per = (100*mean(se_surv)/10)) %>%
  mutate(upper = mean_precip_per + mean_se_per,
         lower = mean_precip_per - mean_se_per) %>% 
  ggplot(mapping = aes(x = precip, y = mean_precip_per, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30","red1", "blue1")) +
  scale_x_discrete(labels = c("Ambient", "Drought", "Wet")) +
  labs(y = "Mean Survival (%)",
       x = "Precipitation Treatment") +
  labs_pubr() +
  theme_pubr(legend = "right")

# plot survival over time with error bars

labels.precip <- factor(tot_surv$precip, labels = c("Ambient", "Drought", "Wet"))

tot_surv %>% 
  mutate(precip = recode(precip,
    "Control" = "Ambient",
    "RO" = "Drought",
    "IR" = "Wet")) %>% 
  ggplot(aes(x = date, y = 10*mean_surv, group = cohort, color = precip)) + 
  scale_color_manual(values = c("grey30", "red1", "blue1")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_line(aes(linetype = cohort), stat = "identity", size = 1) + 
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  ylim(0,100) +
  labs(y = "Mean Survival (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       linetype = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip, ncol = 3, nrow = 1) +
  theme_pubr(legend = "bottom", x.text.angle = 45)

# create data set for precip and excl
tot_surv_pe <- seedlings_obs %>% 
  group_by(precip, excl) %>% 
  summarise(mean_surv = mean(survival),
            se_surv = (sd(survival)/n())) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

# if want to look at effects treating precip as continuous
precip_cont_df <- seedlings_obs %>% 
  ungroup() %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "100",
                                     "IR" = "165",
                                     "RO" = "35")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

# clip vs unclip across exclusion treatments
precip_cont_df %>%
  ggplot(aes(x = precip_cont, y = survival, color = clip, group = clip)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))+
  facet_wrap(~excl)

# only exclusion txs
precip_cont_df %>% 
  ggplot(aes(x = precip_cont, y = survival, color = excl, group = excl)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))

# year 2 zeros in total exclusion tx really bring down the averages
test %>% filter(cohort == "2") %>% 
  ggplot()+
  geom_histogram(aes(x = survival))+
  facet_wrap(~excl)

### Herbivory visualization
seedlings <- vroom("Data/seedlings_combined.csv",
                   col_select = -c(1),
                   col_types = c(.default = "f",
                                 date = "D"))
str(seedlings)

# Herbivory counts
seedlings_herb_full <- seedlings %>% 
  group_by(block, precip, clip, excl, side, rep, date, cohort) %>%
  count(herbivory) %>% 
  pivot_wider(names_from = herbivory,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_herb = "0",
         herb_lived = "1",
         herb_died = "2") %>% 
  mutate(tot_herbivory = herb_lived + herb_died)

glimpse(seedlings_herb_full)
summary(seedlings_herb_full)

# descriptive stats for each cohort 1-3
describeBy(seedlings_herb_full, group = "cohort")

# make a new ID column for each observation and plot
seedlings_obs_herb <- seedlings_herb_full %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("sampID", block:rep, sep = "_", remove = FALSE) %>%
  unite("plotID", block:precip, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         date = as.factor(date),
         sampID = as.factor(sampID),
         plotID = as.factor(plotID))

seedlings_obs_herb %>% 
  mutate(date = as.Date(date)) %>% 
  #select(c(precip, clip, excl, cohort, date, tot_herbivory)) %>% 
  ggplot(aes(x = date, y=herb_died,group = cohort, color = precip)) +
  geom_point(size=2, alpha=0.4) +
  geom_smooth(method="loess", colour="blue", size=1.5) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  scale_color_manual(values = c("grey30", "blue1", "red1")) +
  ylim(0,1)+
  xlab("Date") +
  ylab("Tot_herb") +
  facet_wrap(~precip + excl, scales = "free_y") +
  theme_pubr(legend = "bottom", x.text.angle = 45)

# use summary data
seedlings_obs %>%
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(precip = recode_factor(precip,
                         "Control" = "Ambient",
                         "RO" = "Drought",
                         "IR" = "Wet")) %>%
  ggplot(aes(x = date, y = 100*(tot_herbivory/10), group = cohort, color = precip, linetype = cohort)) + 
  scale_color_manual(values = c("grey30", "red1", "blue1")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  #geom_point() +
  stat_smooth(method = "loess") +
  #geom_line(aes(linetype = cohort), stat = "identity", size = 1) + 
  #scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  #scale_fill_manual(values = c("grey30","blue1", "red1")) +
  #ylim(0, 20) +
  labs(y = "Total Herbivory (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       points = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip + excl, scales = "free_y") +
  theme_pubr(legend = "bottom", x.text.angle = 45)


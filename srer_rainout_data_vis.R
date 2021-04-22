# USDA Rainout - Santa Rita Experimental Range
# Data visualizations
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-03-29

# Load packages
library(tidyverse)
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
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  scale_x_discrete(labels = c("Ambient", "Drought", "Wet")) +
  labs(y = "Mean Survival (%)",
       x = "Precipitation Treatment") +
  labs_pubr() +
  theme_pubr(legend = "right")

# plot predicted survival over time with error bars

labels.precip <- factor(tot_surv$precip, labels = c("Ambient", "Drought", "Wet"))

tot_surv %>% 
  mutate(precip = recode(precip,
    "Control" = "Ambient",
    "RO" = "Drought",
    "IR" = "Wet")) %>% 
  ggplot(aes(x = date, y = 10*mean_surv, group = cohort, color = precip)) + 
  scale_color_manual(values = c("grey30", "red1", "blue1")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_line(aes(linetype = cohort), stat = "identity", size = 2) + 
  #geom_errorbar(aes(ymin = 10*lower, ymax = 10*upper), width = 1, size = 5) + # very busy
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  labs(y = "Mean Survival (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       linetype = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip, ncol = 1, nrow = 3) +
  theme_pubr(legend = "bottom")

# create data set for precip and excl
tot_surv_pe <- seedlings_obs %>% 
  group_by(precip, excl) %>% 
  summarise(mean_surv = mean(survival),
            se_surv = (sd(survival)/n())) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

test <- seedlings_obs %>% 
  ungroup() %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "100",
                                     "IR" = "165",
                                     "RO" = "35")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

test %>%
  ggplot(aes(x = precip_cont, y = survival, color = clip, group = clip)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))+
  facet_wrap(~excl)

test %>% 
  ggplot(aes(x = precip_cont, y = survival, color = excl, group = excl)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))

# year 2 zeros in total exclusion tx really bring down the averages
test %>% filter(cohort == "2") %>% 
  ggplot()+
  geom_histogram(aes(x = survival))+
  facet_wrap(~excl)

# plot control raw survival values for cohort 1
surv_1 <- seedlings_obs %>%
  filter(cohort == '1'& precip == 'Control') %>%
  group_by(date) %>% 
  summarise(mean_surv_1 = mean(survival),
            se_surv = (sd(survival)/n())) %>%
  mutate(upper = mean_surv_1 + se_surv,
         lower = mean_surv_1 - se_surv)

surv_1 %>% 
  ggplot(aes(x = date, y = mean_surv_1)) +
  geom_line(group = 1, color = "grey30")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1)+
  theme_pubr()

# plot control raw survival values for cohort 2
surv_2 <- seedlings_obs %>%
  filter(cohort == '2'& precip == 'Control') %>%
  group_by(date) %>% 
  summarise(mean_surv_2 = mean(survival),
            se_surv = (sd(survival)/n())) %>%
  mutate(upper = mean_surv_2 + se_surv,
         lower = mean_surv_2 - se_surv)

surv_2 %>% 
  ggplot(aes(x = date, y = mean_surv_2)) +
  geom_line(group = 1, color = "grey30")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1)+
  theme_pubr()

# plot control predicted survival values for cohort 3
surv_3 <- seedlings_obs %>%
  filter(cohort == '3'& precip == 'Control') %>%
  group_by(date) %>% 
  summarise(mean_surv_3 = mean(survival),
            se_surv = (sd(survival)/n())) %>%
  mutate(upper = mean_surv_3 + se_surv,
         lower = mean_surv_3 - se_surv)

surv_3 %>% 
  ggplot(aes(x = date, y = mean_surv_3)) +
  geom_line(group = 1, color = "grey30")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1)+
  theme_pubr()

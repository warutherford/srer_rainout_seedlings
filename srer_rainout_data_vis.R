# USDA Rainout - Santa Rita Experimental Range
# Data visualizations-main
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-03-29

# Load packages
library(tidyverse)
library(forcats)
library(vroom)
library(ggpubr)
library(glmmTMB)
library(dplyr)

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

bar_surv_fig <- tot_surv %>%
  group_by(precip) %>% 
  summarise(mean_precip_per = (100*mean(mean_surv)/10),
            mean_se_per = (100*mean(se_surv)/10)) %>%
  mutate(upper = mean_precip_per + mean_se_per,
         lower = mean_precip_per - mean_se_per) %>%
  mutate(precip = recode_factor(precip,
                                       "IR" = "Wet",
                                      "Control" = "Ambient",
                                       "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = mean_precip_per, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("blue1","grey30","#ba7525")) +
  scale_x_discrete(labels = c("Wet","Ambient", "Drought")) +
  labs(y = "Mean Survival (%)",
       x = "Precipitation Treatment") +
  theme_pubr(legend = "none") +
  labs_pubr()

bar_surv_fig

ggsave(filename = "Figures_Tables/bar_mean_surv.tiff",
       plot = bar_surv_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# plot survival over time with error bars

#labels.precip <- factor(tot_surv$precip, labels = c("Ambient", "Drought", "Wet"))

line_mean_surv <- tot_surv %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(aes(x = date, y = 10*mean_surv, group = cohort, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_y_continuous(breaks = seq(0, 110, 10), expand = c(0.02,0)) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_line(aes(linetype = cohort), stat = "identity", size = 1) + 
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  labs(y = "Mean Survival (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       linetype = "Cohort") +
  facet_wrap(~precip, ncol = 3, nrow = 1) +
  theme_pubr(legend = "bottom", x.text.angle = 45) +
  labs_pubr()

line_mean_surv

ggsave(filename = "Figures_Tables/seedlings/line_mean_surv.tiff",
       plot = line_mean_surv,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

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
  group_by(precip_cont, precip, clip, excl) %>% 
  summarise(mean_surv = 100*mean(survival/10)) %>% #convert survival to a percentage
  ggplot(aes(x = precip_cont, y = mean_surv, color = clip, group = clip)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))+
  labs(y = "Mean Survival (%)",
       x = "Precipitation (mm)",
       color = "Grazing") +
  labs_pubr() +
  facet_wrap(~excl, nrow = 1)

# only exclusion txs
precip_cont_df %>% 
  ggplot(aes(x = precip_cont, y = survival, color = excl, group = excl)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))

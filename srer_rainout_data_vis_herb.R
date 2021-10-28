# USDA Rainout - Santa Rita Experimental Range
# Data visualizations-herbivory
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-04-28

# Load packages
library(tidyverse)
library(forcats)
library(vroom)
library(ggpubr)
library(glmmTMB)
library(dplyr)

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

###
## Figures ##
###

# Presence/Absence (0/1) of Herbivory that resulting in seedling death across all years
died_herb <- seedlings_obs_herb %>% 
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(date = as.Date(date)) %>% 
  ggplot(aes(x = date, y = herb_died, group = precip, color = precip)) +
  geom_point(size = 2, alpha = 0.4) +
  geom_smooth(method = "loess", color = "#885a99", size = 1.5, se = FALSE) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  ylim(0, 1) +
  xlab("Date") +
  ylab("Died Following Herbivory") +
  facet_wrap(~precip, ncol = 3, scales = "free_y") +
  theme_pubr(legend = "bottom", x.text.angle = 45)

died_herb

# Use summary data
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

# Percentage of seedlings that died from herbivory by pptx, cohort/year, and clipping
herb_died_time_clip_fig <- seedlings_obs %>%
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "RO" = "Drought",
                                "IR" = "Wet")) %>%
  ggplot(aes(x = date, y = 100*(herb_died/10), group = cohort, color = precip, linetype = cohort)) + 
  scale_color_manual(values = c("grey30", "#ba7525", "blue1")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  stat_smooth(method = "loess") +
  #ylim(0, 15) +
  labs(y = "Died Following Herbivory (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       points = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip + clip, ncol = 2) +
  theme_pubr(legend = "bottom", x.text.angle = 45)

ggsave(filename = "Figures_Tables/herb_died_ppt_clip.tiff",
       plot = herb_died_time_clip_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# Percentage of seedlings that died from herbivory by pptx, cohort/year
herb_died_time_fig <- seedlings_obs %>%
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "RO" = "Drought",
                                "IR" = "Wet")) %>%
  ggplot(aes(x = date, y = 100*(herb_died/10), group = cohort, color = precip, linetype = cohort)) + 
  scale_color_manual(values = c("grey30", "#ba7525", "blue1")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  stat_smooth(method = "loess") +
  #ylim(0, 15) +
  labs(y = "Died Following Herbivory (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       points = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip+clip, ncol = 2) +
  theme_pubr(legend = "bottom", x.text.angle = 45)

ggsave(filename = "Figures_Tables/herb_died_ppt.tiff",
       plot = herb_died_time_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# create herbivory only df for quicker graphing
herb_df_all <- seedlings_obs %>% 
  group_by(precip, clip, excl) %>% 
  summarise(died_mean = 100*(mean(herb_died)/10),
            live_mean = 100*(mean(herb_lived)/10),
            tot_mean = 100*(mean(tot_herbivory)/10),
            mean_herbD_se = 100*(sd(herb_died)/sqrt(n())),
            mean_herbL_se = 100*(sd(herb_lived)/sqrt(n())),
            mean_herbTOT_se = 100*(sd(tot_herbivory)/sqrt(n()))) %>% 
  mutate(upper_herbD = died_mean + mean_herbD_se,
         lower_herbD = died_mean - mean_herbD_se,
         upper_herbL = live_mean + mean_herbL_se,
         lower_herbL = live_mean - mean_herbL_se,
         upper_herbTOT = tot_mean + mean_herbTOT_se,
         lower_herbTOT = tot_mean - mean_herbTOT_se) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Rodents Excl",
                              "Total" = "Total Excl")) %>% 
  ungroup()

# bar graph of total percentage of seedlings with herbivory
herb_total_fig <- herb_df_all %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = tot_mean, fill = precip)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower_herbTOT, ymax = upper_herbTOT), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 35) +
  labs(y = "Total Herbivory (%)",
       x = "PPTx") +
  labs_pubr() +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl))

ggsave(filename = "Figures_Tables/herb_total_bar.tiff",
       plot = herb_total_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# bar graph of percentage of seedlings that died from herbivory
herb_died_fig <- herb_df_all %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = died_mean, fill = precip)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower_herbD, ymax = upper_herbD), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 30) +
  labs(y = "Died Following Herbivory (%)",
       x = "PPTx") +
  labs_pubr() +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl))

ggsave(filename = "Figures_Tables/herb_died_bar.tiff",
       plot = herb_died_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# bar graph of percentage of seedlings that lived following herbivory
herb_lived_fig <- herb_df_all %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = live_mean, fill = precip)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower_herbL, ymax = upper_herbL), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  labs(y = "Lived Following Herbivory (%)",
       x = "PPTx") +
  #ylim(0, 5)+
  labs_pubr() +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl))

ggsave(filename = "Figures_Tables/herb_lived_bar.tiff",
       plot = herb_lived_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# use change from previous measurement of herbivory presence/absence (0/1)
change <- seedlings_obs_herb %>% 
  group_by(precip, clip, excl, side, rep, date) %>% 
  mutate(date = as.Date(date)) %>% 
  arrange(desc(date), .by_group = TRUE) %>% 
  mutate(herb_died_lag = dplyr::lag(herb_died, n = 1, default = NA, order_by = date),
         herb_lived_lag = dplyr::lag(herb_lived, n = 1, default = NA, order_by = date),
         herb_tot_lag = dplyr::lag(tot_herbivory, n = 1, default = NA, order_by = date))

# below plots all herbivory regardless of resulting death or live across precip and clipping txs
died_herb_change_fig <- change %>% 
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought"),
         date = as.Date(date)) %>% 
  ggplot(aes(x = date, y = 100*herb_died_lag, color = precip, group = cohort, linetype = cohort)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "4 months") +
  stat_smooth(method = "loess", span = 0.25) +
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  labs(y = "Change in Seedling Death Following Herbivory (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       linetype = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip + clip, ncol = 2) +
  theme_pubr(legend = "bottom", x.text.angle = 45, border = TRUE)
# 16080 values are NAs (missing) because of the lag calculation

ggsave(filename = "Figures_Tables/change_herb_died.tiff",
       plot = died_herb_change_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# below plots all herbivory regardless of resulting in death or live across precip and clipping txs
tot_herb_change_fig <- change %>% 
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought")) %>% 
  ggplot(aes(x = date, y = 100*herb_tot_lag, color = precip, group = cohort, linetype = cohort)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  stat_smooth(method = "loess", span = 0.25) +
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  labs(y = "Change in Seedling Herbivory (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       linetype = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip + clip, ncol = 2) +
  theme_pubr(legend = "bottom", x.text.angle = 45)

ggsave(filename = "Figures_Tables/change_herb_total.tiff",
       plot = tot_herb_change_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# below plots all herbivory regardless of resulting death or live, but only
# the patterns of the precip treatments (no clipping or exclusion)
# figure takes a lot of time to run!
tot_herb_change_fig_all <- change %>% 
  group_by(precip, date) %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought")) %>% 
  ggplot(aes(x = date, y = 100*herb_tot_lag, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_smooth(method = "loess", se = TRUE, span = 0.25) +
  ylim(0,NA) +
  labs(y = "Change in Seedling Herbivory (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       points = "Cohort") +
  labs_pubr() +
  facet_wrap(~precip, ncol = 3) +
  theme_pubr(legend = "bottom", x.text.angle = 45, border = TRUE)

ggsave(filename = "Figures_Tables/change_herb_total_ppt.tiff",
       plot = tot_herb_change_fig_all,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


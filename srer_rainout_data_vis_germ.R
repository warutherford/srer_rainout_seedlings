# USDA Rainout - Santa Rita Experimental Range
# Data visualizations-germination
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-05-05

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

# For germination data, either the seed germinated or it didn't
seedlings_germ_full <- seedlings %>% 
  group_by(block, precip, clip, excl, side, rep, date, cohort) %>%
  count(fate) %>% 
  pivot_wider(names_from = fate,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_germ = "0",
         survival = "1",
         died = "2") %>% 
  mutate(tot_germination = (survival + died))

# Create possible random factor variables and fix any data structures needed for modeling
seedlings_obs_germ <- seedlings_germ_full %>% 
  mutate(block = as.character(block),
         precip = as.character(precip)) %>% 
  unite("plotID",block:precip, sep = "_", remove = FALSE) %>%
  unite("sampID", block:rep, sep = "_", remove = FALSE) %>% 
  mutate(block = as.factor(block),
         precip = as.factor(precip),
         plotID = as.factor(plotID),
         date = as.factor(date),
         sampID = as.factor(sampID))

# descriptive stats for each cohort 1-3
describeBy(seedlings_obs_germ , group = "cohort")

# mixed effects model with nesting, sampID as random and date (cohort) within year for temporal autocorrelation
hist(seedlings_obs_germ$tot_germination)# binomial (0 or 1)

###
## Figures ##
###

# Presence/Absence (0/1) of germination
tot_germ_fig <- seedlings_obs_germ %>% 
  group_by(precip, clip, excl, cohort, date) %>% 
  mutate(date = as.Date(date)) %>% 
  ggplot(aes(x = date, y = tot_germination, group = precip, color = precip)) +
  geom_point(size=2, alpha=0.4) +
  geom_smooth(method="loess", colour="blue", size=1.5, se = FALSE) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  ylim(0,1)+
  xlab("Date") +
  ylab("Germination") +
  facet_wrap(~precip, ncol = 3, scales = "free_y") +
  theme_pubr(legend = "bottom", x.text.angle = 45) +
  labs_pubr()

tot_germ_fig

# Use summary data
### Read in clean seedlings data with model outputs from "srer_rainout_stats_summarize"
# seedlings_obs <- vroom("Data/seedlings_obs.csv",
#                        col_types = c(.default = "f",
#                                      ObsID = "i",
#                                      date = "D",
#                                      survival = "i",
#                                      died = "i",
#                                      tot_germination = "i",
#                                      herb_lived = "i",
#                                      herb_died = "i",
#                                      tot_herbivory = "i",
#                                      granivory = "i",
#                                      res_surv = "d",
#                                      sim_surv = "d",
#                                      sim_fit_surv = "d",
#                                      pred_surv = "d"))
# str(seedlings_obs)
# glimpse(seedlings_obs)


# create data set for precip and excl and clip
tot_germ_pce <- seedlings_obs %>% 
  group_by(precip, excl, clip) %>% 
  summarise(mean_germ = 100*mean(tot_germination/10),
            sd_germ = 100*sd(tot_germination/10),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + se_germ,
         lower = mean_germ - se_germ)

bar_pce_germ_fig <- tot_germ_pce %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Rodents Excl",
                              "Total" = "Total Excl")) %>% 
  ggplot(mapping = aes(x = precip, y = mean_germ, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 80) +
  labs(y = "Mean Germination (%)",
       x = "PPTx") +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl)) +
  labs_pubr()

bar_pce_germ_fig

ggsave(filename = "Figures_Tables/bar_alltx_germ.tiff",
       plot = bar_pce_germ_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# create data set for precip and clip
tot_germ_pc <- seedlings_obs %>% 
  group_by(precip, clip, cohort) %>% 
  summarise(mean_germ = 100*mean(tot_germination/10),
            sd_germ = 100*sd(tot_germination/10),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + se_germ,
         lower = mean_germ - se_germ)

bar_clip_germ_fig <- tot_germ_pc %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = clip, y = mean_germ, fill = clip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("brown","darkorange")) +
  scale_x_discrete(labels = c("Clipped","Unclipped")) +
  ylim(0, 80) +
  labs(y = "Mean Germination (%)",
       x = "Grazing Treatment") +
  theme_pubr(legend = "none") +
  facet_wrap(~precip+cohort)+
  labs_pubr()

bar_clip_germ_fig

ggsave(filename = "Figures_Tables/bar_clip_germ.tiff",
       plot = bar_clip_germ_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# create data set for clip
tot_surv_c <- seedlings_obs %>% 
  group_by(clip,cohort) %>% 
  summarise(mean_germ = 100*mean(tot_germination/10),
            sd_germ = 100*sd(tot_germination/10),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + se_germ,
         lower = mean_germ - se_germ)

bar_cliponly_germ_fig <- tot_surv_c %>%
  ggplot(mapping = aes(x = clip, y = mean_germ, fill = clip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("brown","darkorange")) +
  scale_x_discrete(labels = c("Clipped","Unclipped")) +
  labs(y = "Mean Germination (%)",
       x = "Grazing Treatment") +
  theme_pubr(legend = "none") +
  facet_wrap(~cohort)+
  labs_pubr()

bar_cliponly_germ_fig

ggsave(filename = "Figures_Tables/bar_clipcohort_germ.tiff",
       plot = bar_cliponly_germ_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

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
  summarise(mean_germ = 100*mean(tot_germination/10)) %>% #convert survival to a percentage
  ggplot(aes(x = precip_cont, y = mean_germ, color = clip, group = clip)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))+
  labs(y = "Mean Germination (%)",
       x = "Precipitation (mm)",
       color = "Grazing") +
  labs_pubr() +
  facet_wrap(~excl, nrow = 1)

# only exclusion txs
precip_cont_df %>% 
  ggplot(aes(x = precip_cont, y = 10*tot_germination, color = excl, group = excl)) + 
  geom_smooth(method = "glm", formula = y ~ log(x))


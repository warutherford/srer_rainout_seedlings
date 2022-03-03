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

seedlings_obs_germ <- vroom("Data/seedlings_obs_germ.csv",
                   col_select = -c(1),
                   col_types = c(.default = "d",
                                 sampID = "f",
                                 block = "f",
                                 date = "D",
                                 precip = "f",
                                 clip = "f",
                                 excl = "f",
                                 side = "f",
                                 rep = "f",
                                 cohort = "f"))
str(seedlings_obs_germ)

# descriptive stats for each cohort 1-3
describeBy(seedlings_obs_germ , group = "cohort")

# mixed effects model with nesting, sampID as random and date (cohort) within year for temporal autocorrelation
hist(seedlings_obs_germ$tot_germination)# binomial (0 or 1)

# summary for germ by only precip
tot_germ_p <- seedlings_obs_germ %>% 
  group_by(precip) %>% 
  summarise(mean_germ = 100*mean(tot_germination),
            sd_germ = 100*sd(tot_germination),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + se_germ,
         lower = mean_germ - se_germ)

tot_germ_p

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


# create data set for precip and excl and clip
tot_germ_pce <- seedlings_obs_germ %>% 
  group_by(precip, excl, clip) %>% 
  summarise(mean_germ = mean(100*(tot_germination)),
            sd_germ = sd(100*(tot_germination)),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + 10*se_germ,
         lower = mean_germ - 10*se_germ)

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
  labs(y = "Seed Germination (%)",
       x = "Precipitation Treatment") +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl)) +
  labs_pubr(base_size = 24) +
  theme(legend.position="none", 
        panel.border = element_blank(), 
        panel.spacing.x = unit(1,"line"))

bar_pce_germ_fig

ggsave(filename = "Figures_Tables/bar_alltx_germ.tiff",
       plot = bar_pce_germ_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


# create data set for precip and clip
tot_germ_pc <- seedlings_obs_germ %>% 
  group_by(precip, clip) %>% 
  summarise(mean_germ = 100*mean(tot_germination),
            sd_germ = 100*sd(tot_germination),
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
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10)) +
  labs(y = "Mean Germination (%)",
       x = "") +
  theme_pubr(legend = "none") +
  facet_wrap(~precip)+
  labs_pubr(base_size = 24)

bar_clip_germ_fig

ggsave(filename = "Figures_Tables/bar_clip_germ.tiff",
       plot = bar_clip_germ_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# create data set for clip
tot_germ_c <- seedlings_obs_germ %>% 
  group_by(clip,cohort) %>% 
  summarise(mean_germ = 100*mean(tot_germination),
            sd_germ = 100*sd(tot_germination),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + se_germ,
         lower = mean_germ - se_germ)

bar_cliponly_germ_fig <- tot_germ_c %>%
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



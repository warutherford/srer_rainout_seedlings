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

# create data set for precip and excl and cohort
tot_germ_pce <- seedlings_obs_germ %>% 
  group_by(precip, clip, excl, cohort) %>% 
  summarise(mean_germ = mean(100*(tot_germination)),
            sd_germ = sd(100*(tot_germination)),
            counts = n(),
            se_germ = (sd_germ/sqrt(counts))) %>%
  mutate(upper = mean_germ + 10*se_germ,
         lower = mean_germ - 10*se_germ)

bar_pce_germ_fig_c <- tot_germ_pce %>%
  filter(clip == "Clipped") %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "All Excl",
                              "Ants" = "Rodents Excl",
                              "Rodents" = "Ants Excl",
                              "Total" = "None")) %>% 
  ggplot(mapping = aes(x = precip, y = mean_germ, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 100) +
  labs(y = "Seed Germination (%)",
       x = "Precipitation Treatment") +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(cohort), rows = vars(excl)) +
  labs_pubr(base_size = 24) +
  theme(legend.position="none", 
        panel.border = element_blank(), 
        panel.spacing.x = unit(1,"line"))

bar_pce_germ_fig_c

ggsave(filename = "Figures_Tables/bar_peco_germ_clip_fig2.tiff",
       plot = bar_pce_germ_fig_c,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

bar_pce_germ_fig_uc <- tot_germ_pce %>%
  filter(clip == "Unclipped") %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "All Excl",
                              "Ants" = "RodentsExcl",
                              "Rodents" = "Ants Excl",
                              "Total" = "None")) %>% 
  ggplot(mapping = aes(x = precip, y = mean_germ, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 100) +
  labs(y = "Seed Germination (%)",
       x = "Precipitation Treatment") +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(cohort), rows = vars(excl)) +
  labs_pubr(base_size = 24) +
  theme(legend.position="none", 
        panel.border = element_blank(), 
        panel.spacing.x = unit(1,"line"))

bar_pce_germ_fig_uc

ggsave(filename = "Figures_Tables/bar_peco_germ_unclip_fig2.tiff",
       plot = bar_pce_germ_fig_uc,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")





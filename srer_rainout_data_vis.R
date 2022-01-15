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
            sd_surv = sd(survival),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv,
         date = as.Date(date))

summary(tot_surv)


bar_surv_fig <- tot_surv %>%
  group_by(precip) %>% 
  summarise(mean_precip_per = (100*mean(mean_surv)/10),
            mean_se_per = (100*mean(se_surv)/10)) %>%
  mutate(upper = mean_precip_per + mean_se_per,
         lower = mean_precip_per - mean_se_per) %>%
  mutate(precip = recode_factor(precip,
                                      "Control" = "Ambient",
                                      "IR" = "Wet",
                                       "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = mean_precip_per, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30","blue1","#ba7525")) +
  scale_x_discrete(labels = c("Ambient","Wet", "Drought")) +
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

# create data set for clip
tot_surv_c <- seedlings_obs %>% 
  group_by(clip) %>% 
  summarise(mean_surv = 100*mean(survival/10),
            sd_surv = 100*sd(survival/10),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

bar_cliponly_fig <- tot_surv_c %>%
  ggplot(mapping = aes(x = clip, y = mean_surv, fill = clip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("brown","darkorange")) +
  scale_x_discrete(labels = c("Clipped","Unclipped")) +
  labs(y = "Mean Survival (%)",
       x = "Grazing Treatment") +
  theme_pubr(legend = "none") +
  labs_pubr()

bar_cliponly_fig

ggsave(filename = "Figures_Tables/bar_cliponly_surv.tiff",
       plot = bar_cliponly_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# plot survival over time, all dates
line_mean_surv <- tot_surv %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(aes(x = date, y = 10*mean_surv, group = cohort, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  #scale_y_continuous(breaks = seq(0, 120, 10), expand = c(0.02,0)) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_line(aes(linetype = cohort), stat = "identity", size = 1) + 
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  ylim(0, 100) +
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


# remove the lead up survival to start max for each cohort (following Archer suggestions)
surv_small_1 <- tot_surv %>% 
  filter(cohort == 1) %>% 
  filter(date > "2017-07-21")

surv_small_2 <- tot_surv %>% 
  filter(cohort == 2) %>% 
  filter(date >= "2018-07-18")

surv_small_3 <- tot_surv %>% 
  filter(cohort == 3 & precip != "RO") %>% 
  filter(date >= "2019-08-03")

surv_small_drought <- tot_surv %>% 
  filter(precip == "RO" & date >= "2019-08-10")

surv_small_all <- rbind(surv_small_1, surv_small_2, surv_small_3, surv_small_drought)

line_mean_surv_archer <- surv_small_all %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(aes(x = date, y = 10*mean_surv, group = cohort, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  #scale_y_continuous(breaks = seq(0, 120, 10), expand = c(0.02,0)) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_line(aes(linetype = cohort), stat = "identity", size = 1) + 
  scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  ylim(0, 100) +
  labs(y = "Mean Survival (%)",
       x = "Date (Month-Year)",
       color = "PPTx",
       linetype = "Cohort") +
  facet_wrap(~precip, ncol = 3, nrow = 1) +
  theme_pubr(legend = "bottom", x.text.angle = 45) +
  labs_pubr()

line_mean_surv_archer

ggsave(filename = "Figures_Tables/seedlings/line_mean_surv_archer.tiff",
       plot = line_mean_surv_archer,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# calculate summary info
tot_surv_summary <- seedlings_obs %>% 
  group_by(precip, cohort, date) %>% 
  summarise(mean_surv = 100*mean(survival/10),
            sd_surv = 100*sd(survival/10),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_summary_full <- seedlings_obs %>% 
  group_by(precip, excl, clip) %>% 
  summarise(mean_surv = 100*mean(survival/10),
            sd_surv = 100*sd(survival/10),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_summary_pcohort <- seedlings_obs %>% 
  group_by(precip, cohort) %>% 
  summarise(mean_surv = 100*mean(survival/10),
            sd_surv = 100*sd(survival/10),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

# create data set for precip and excl and clip
tot_surv_pce <- seedlings_obs %>% 
  group_by(precip, excl, clip) %>% 
  summarise(mean_surv = 100*mean(survival/10),
            sd_surv = 100*sd(survival/10),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

bar_pce_fig <- tot_surv_pce %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Rodents Excl",
                              "Total" = "Total Excl")) %>% 
  ggplot(mapping = aes(x = precip, y = mean_surv, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 35) +
  labs(y = "Seedling Survival (%)",
       x = "PPTx") +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl)) +
  labs_pubr(base_size = 24) +
  theme(legend.position="none", 
        panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"))

bar_pce_fig

ggsave(filename = "Figures_Tables/bar_alltx_surv.tiff",
       plot = bar_pce_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# create data set for precip and clip
tot_surv_pc <- seedlings_obs %>% 
  group_by(precip, clip) %>% 
  summarise(mean_surv = 100*mean(survival/10),
            sd_surv = 100*sd(survival/10),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

bar_clip_fig <- tot_surv_pc %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = clip, y = mean_surv, fill = clip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("brown","darkorange")) +
  scale_x_discrete(labels = c("Clipped","Unclipped")) +
  ylim(0, 30) +
  labs(y = "Seedling Survival (%)",
       x = "") +
  theme_pubr(legend = "none") +
  facet_wrap(~precip)+
  labs_pubr(base_size = 24)

bar_clip_fig

ggsave(filename = "Figures_Tables/bar_clip_surv.tiff",
       plot = bar_clip_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# if want to look at effects treating monsoon precip as continuous
precip_cont_surv_df_1 <- seedlings_obs %>% 
  group_by(cohort) %>%
  filter(cohort == "1") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "280",
                                     "IR" = "462",
                                     "RO" = "98")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_surv_df_2 <-seedlings_obs %>% 
  group_by(cohort) %>%
  filter(cohort == "2") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "330",
                                     "IR" = "545",
                                     "RO" = "115")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_surv_df_3 <-seedlings_obs %>% 
  group_by(cohort) %>%
  filter(cohort == "3") %>% 
  mutate(precip_cont = dplyr::recode(precip,
                                     "Control" = "292",
                                     "IR" = "482",
                                     "RO" = "102")) %>% 
  mutate(precip_cont = as.numeric(as.character(precip_cont)))

precip_cont_surv_df <- rbind(precip_cont_surv_df_1, precip_cont_surv_df_2, precip_cont_surv_df_3)


# predicted survival only PPT
pred_surv_pt <- precip_cont_surv_df %>% 
  group_by(precip_cont, clip, excl) %>% 
  summarise(mean_pred_surv = mean(pred_surv),
            sd_surv = sd(pred_surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_pred_surv + se_surv,
         lower = mean_pred_surv - se_surv) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Small Mammals Excl",
                              "Total" = "All Excl")) %>% 
  ggplot(aes(x = precip_cont, y = 10*(mean_pred_surv))) + 
  #geom_point() +
  geom_pointrange(aes(ymin = 10*lower, ymax = 10*upper, color = excl), size = 0.5) +
  geom_smooth(method = "glm", formula = y ~ log(x) + x, se = T, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Precipitation (mm)",
       color = "Exclusion") +
  scale_x_continuous(breaks = c(0,50, 100,150, 200,250, 300,350, 400,450, 500, 550), limits = c(0, 550))+
  ylim(0, 45) +
  theme_pubr(legend = "right")+
  labs_pubr(base_size = 24)

pred_surv_pt

ggsave(filename = "Figures_Tables/pred_surv_cont.tiff",
       plot = pred_surv_pt,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


pred_surv_excl <- precip_cont_surv_df %>% 
  group_by(precip_cont, excl) %>% 
  summarise(mean_pred_surv = mean(pred_surv),
            sd_surv = sd(pred_surv),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_pred_surv + se_surv,
         lower = mean_pred_surv - se_surv) %>% 
  mutate(excl = recode_factor(excl, 
                              "Control" = "None",
                              "Ants" = "Ants Excl",
                              "Rodents" = "Rodents Excl",
                              "Total" = "All Excl")) %>% 
  ggplot(aes(x = precip_cont, y = 10*(mean_pred_surv), group = excl, color = excl)) + 
  geom_smooth(method = "glm", formula = y ~ log(x) + x, se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Precipitation (mm)",
       color = "Exclusion") +
  scale_x_continuous(breaks = c(0,50, 100,150, 200,250, 300,350, 400,450, 500, 550), limits = c(0, 550))+
  ylim(0, 45) +
  theme_pubr(legend = "none")+
  labs_pubr(base_size = 24)

pred_surv_excl

ggsave(filename = "Figures_Tables/pred_surv_cont_excl.tiff",
       plot = pred_surv_excl,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


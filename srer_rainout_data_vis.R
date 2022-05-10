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
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, cohort, date) %>% 
  summarise(mean_surv = mean(survival/10),
            count_surv = sum(survival),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv,
         date = as.Date(date))

summary(tot_surv)

# survivorship curve, remove the lead up
surv_small_1 <- tot_surv %>% 
  filter(cohort == 1) %>% 
  filter(date > "2017-07-21")

surv_small_2 <- tot_surv %>% 
  filter(cohort == 2) %>% 
  filter(date >= "2018-07-18")

surv_small_3 <- tot_surv %>% 
  filter(cohort == 3 & precip != "RO") %>% # to start with all max values, need to cut one more date from RO..delayed germ
  filter(date >= "2019-08-03")

surv_small_drought <- tot_surv %>% 
  filter(precip == "RO" & date >= "2019-08-10")

surv_small_all <- rbind(surv_small_1, surv_small_2, surv_small_3, surv_small_drought) %>% 
  mutate(date = as.Date(date))


line_surv_curv <- surv_small_all %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought")) %>% 
  ggplot(aes(x = date, y = mean_surv, group = precip, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "3 months") +
  geom_line(aes(), stat = "identity", size = 1) + 
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  scale_y_log10(breaks = seq(0, 1, 0.10))+
  labs(y = "Proportion Survival",
       x = "Date (Month-Year)",
       color = "PPTx") +
  facet_wrap(~cohort, ncol = 3, nrow = 1) +
  theme_pubr(legend = "bottom", x.text.angle = 45) +
  labs_pubr()

line_surv_curv

# bar fig of survival
bar_surv_fig <- tot_surv %>%
  group_by(precip) %>% 
  summarise(mean_precip_per = (100*mean(mean_surv)),
            mean_se_per = (100*mean(se_surv))) %>%
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
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(clip) %>% 
  summarise(mean_surv = mean(surv_perc),
            count_surv = sum(surv_perc),
            sd_surv = sd(surv_perc),
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
  ggplot(aes(x = date, y = 100*mean_surv, group = cohort, color = precip)) + 
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
  filter(cohort == 3 & precip != "RO") %>% # to start with all max values, need to cut one more date from RO..delayed germ
  filter(date >= "2019-08-03")

surv_small_drought <- tot_surv %>% 
  filter(precip == "RO" & date >= "2019-08-10")

surv_small_all <- rbind(surv_small_1, surv_small_2, surv_small_3, surv_small_drought) %>% 
  mutate(date = as.Date(date))

range <-  c(as.Date("2017-06-15"), as.Date("2020-02-01"))

line_mean_surv_archer <- surv_small_all %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>%
  mutate(cohort = recode_factor(cohort,
                                "1" = "2017 Cohort",
                                "2" = "2018 Cohort",
                                "3" = "2019 Cohort")) %>%
  ggplot(aes(x = date, y = 100*mean_surv, group = precip, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 25), expand = c(0,0)) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "2.5 months", expand = c(0.03,0.03), limits = range) +
  geom_line(aes(), stat = "identity", size = 2.5, position = "jitter") + 
  #scale_linetype_manual(values=c("solid","longdash", "dotted")) +
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  #ylim(0, 700) +
  labs(y = "Survival (%)",
       x = "Date (Month-Year)",
       color = "PPTx") +
  facet_wrap(~cohort, ncol = 3, nrow = 1, scales = "free_x") +
  theme_pubr(legend = "right", x.text.angle = 45) +
  theme(panel.spacing.x = unit(4, "lines")) +
  labs_pubr(base_size = 24)+
  theme(axis.text.x=element_text(size=16))

line_mean_surv_archer

ggsave(filename = "Figures_Tables/seedlings/line_mean_surv_archer.tiff",
       plot = line_mean_surv_archer,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# create a separate plot for each cohort, then combine to set date ranges/scales
# for each cohort 

surv_1 <- surv_small_all %>% filter(cohort == "1")

surv_2 <- surv_small_all %>% filter(cohort == "2")

surv_3 <- surv_small_all %>% filter(cohort == "3")

range1 <-  c(as.Date("2017-05-01"), as.Date("2020-02-01"))

line_mean_surv_archer_1 <- surv_1 %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>%
  mutate(cohort = recode_factor(cohort,
                                "1" = "2017 Cohort")) %>%
  ggplot(aes(x = date, y = 100*mean_surv, group = precip, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 25), expand = c(0.1,0.1)) +
  scale_x_date(date_labels = "%b-%y", date_breaks = "2.75 months", expand = c(0.01, 0), limits = range1) +
  geom_line(aes(), stat = "identity", size = 2.5, position = "jitter") + 
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  labs(y = "Survival (%)",
       x = NULL,
       color = "PPTx") +
  facet_wrap(~cohort) +
  theme_pubr(legend = "none", x.text.angle = 45) +
  #theme(plot.margin =  unit(c(0,1,0,0),"inches")) +
  labs_pubr(base_size = 24)+
  theme(axis.text.x=element_text(size=16))

line_mean_surv_archer_1

range2 <-  c(as.Date("2018-05-01"), as.Date("2020-02-01"))

line_mean_surv_archer_2 <- surv_2 %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>%
  mutate(cohort = recode_factor(cohort,
                                "2" = "2018 Cohort")) %>%
  ggplot(aes(x = date, y = 100*mean_surv, group = precip, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_y_continuous(limits = c(-1,100), breaks = seq(0, 100, 25), expand = c(0.1,0.1)) +
  scale_x_date(date_labels = "%b-%y", date_breaks = "2.75 months", expand = c(0.01,0), limits = range2) +
  geom_line(aes(), stat = "identity", size = 2.5, position = "jitter") + 
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  labs(y = NULL,
       x = NULL,
       color = "PPTx") +
  facet_wrap(~cohort) +
  theme_pubr(legend = "none", x.text.angle = 45) +
  #theme(plot.margin =  unit(c(0,1,0,0),"inches")) +
  labs_pubr(base_size = 24)+
  theme(axis.text.x=element_text(size=16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

line_mean_surv_archer_2

range3 <-  c(as.Date("2019-05-01"), as.Date("2020-02-01"))

line_mean_surv_archer_3 <- surv_3 %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>%
  mutate(cohort = recode_factor(cohort,
                                "3" = "2019 Cohort")) %>%
  ggplot(aes(x = date, y = 100*mean_surv, group = precip, color = precip)) + 
  scale_color_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 25), expand = c(0.1,0.1)) +
  scale_x_date(date_labels = "%b-%y", date_breaks = "2.75 months", expand = c(0.01,0), limits = range3) +
  geom_line(aes(), stat = "identity", size = 2.5, position = "jitter") + 
  scale_fill_manual(values = c("grey30","blue1", "red1")) +
  labs(y = NULL,
       x = NULL,
       color = "PPTx") +
  facet_wrap(~cohort) +
  theme_pubr(legend = "none", x.text.angle = 45) +
  #theme(plot.margin =  unit(c(0,1,0,0),"inches")) +
  labs_pubr(base_size = 24)+
  theme(axis.text.x=element_text(size=16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

line_mean_surv_archer_3

# combine the three figs into one
tot_surv_fig <- gridExtra::grid.arrange(line_mean_surv_archer_1, 
                             line_mean_surv_archer_2,
                             line_mean_surv_archer_3,
                             ncol = 3)

tot_surv_fig

ggsave(filename = "Figures_Tables/seedlings/line_mean_surv_tot.tiff",
       plot = tot_surv_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

surv_small_all %>% 
  group_by(precip, cohort) %>% 
  summarise(max = max(count_surv),
            min = min(count_surv))


# calculate summary info
# remove 2 weeks following sowing
surv_leadup_1 <- seedlings_obs %>% 
  filter(cohort == 1) %>% 
  mutate(date = as.Date(date)) %>% 
  filter(date > "2017-07-21")

surv_leadup_2 <- seedlings_obs %>% 
  filter(cohort == 2) %>% 
  mutate(date = as.Date(date)) %>% 
  filter(date >= "2018-07-18")

surv_leadup_3 <- seedlings_obs%>% 
  filter(cohort == 3 & precip != "RO") %>% # to start with all max values, need to cut one more date from RO..delayed germ
  mutate(date = as.Date(date)) %>% 
  filter(date >= "2019-08-03")

surv_leadup_drought <- seedlings_obs %>% 
  mutate(date = as.Date(date)) %>% 
  filter(precip == "RO" & date >= "2019-08-10")

surv_leadup_all <- rbind(surv_leadup_1, surv_leadup_2, surv_leadup_3, surv_leadup_drought) %>% 
  mutate(date = as.Date(date))


tot_surv_summary <- surv_leadup_all %>% 
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, cohort, date) %>% 
  summarise(mean_surv = 100*mean(surv_perc),
            sd_surv = 100*sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_summary_full <- surv_leadup_all %>% 
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, excl, clip) %>% 
  summarise(mean_surv = 100*mean(surv_perc),
            sd_surv = 100*sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_summary_pcohort <- surv_leadup_all %>% 
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, cohort) %>% 
  summarise(mean_surv = 100*mean(surv_perc),
            sd_surv = 100*sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

# create data set for precip and excl and clip
tot_surv_pce <- surv_leadup_all %>% 
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, excl, clip) %>% 
  summarise(mean_surv = 100*mean(surv_perc),
            sd_surv = 100*sd(surv_perc),
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
  ylim(0, 60) +
  labs(y = "Seedling Survival (%)",
       x = "Precipitation Treatment") +
  theme_pubr(legend = "none") +
  facet_grid(cols = vars(clip), rows = vars(excl)) +
  labs_pubr(base_size = 24) +
  theme(legend.position="none", 
        panel.border = element_blank(), 
        panel.spacing.x = unit(1,"line"))

bar_pce_fig

ggsave(filename = "Figures_Tables/bar_alltx_surv.tiff",
       plot = bar_pce_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# create data set for precip and excl, by cohort
tot_surv_pe <- seedlings_obs %>%
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, excl, cohort) %>% 
  summarise(mean_surv = 100*mean(surv_perc),
            sd_surv = 100*sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

bar_excl_fig <- tot_surv_pe %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = mean_surv, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  #scale_fill_manual(values = c("brown","darkorange")) +
  #scale_x_discrete(labels = c("Clipped","Unclipped")) +
  #ylim(0, 50) +
  labs(y = "Seedling Survival (%)",
       x = "") +
  theme_pubr(legend = "none") +
  facet_grid(rows = vars(cohort), cols = vars(clip), scales = "free") +
  labs_pubr(base_size = 24)

bar_excl_fig


# create data set for precip and clip, by survival year
tot_surv_pc <- seedlings_obs_year %>% # seedslings_obs_year created in 'srer_rainour_stas_summarize' script
  mutate(surv_perc = ((survival/tot_germination))) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0)) %>% 
  group_by(precip, clip, year) %>% 
  summarise(mean_surv = 100*mean(surv_perc),
            sd_surv = 100*sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

bar_clip_fig <- tot_surv_pc %>%
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = precip, y = mean_surv, fill = precip)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(), size = 1) +
  scale_fill_manual(values = c("grey30", "blue1", "#ba7525")) +
  scale_x_discrete(labels = c("Ambient", "Wet", "Drought")) +
  ylim(0, 50) +
  labs(y = "Seedling Survival (%)",
       x = "") +
  theme_pubr(legend = "none") +
  facet_grid(rows = vars(year), cols = vars(clip), scales = "free") +
  labs_pubr(base_size = 24)

bar_clip_fig


ggsave(filename = "Figures_Tables/bar_clip_surv.tiff",
       plot = bar_clip_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")


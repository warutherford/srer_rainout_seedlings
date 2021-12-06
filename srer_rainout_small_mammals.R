# USDA Rainout - Santa Rita Experimental Range
# Small Mammal Trapping Tables
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-12-10

# load packages
library(tidyverse)
library(vroom)
library(gt)
library(glue)

# read in data
rodents <- vroom("Data/Trapping_Combined_2017_2019.csv",
                   col_types = c(.default = "f",
                                 date = "D"))

glimpse(rodents)

rod_new <- rodents %>%
  dplyr::select(-c(date,trap_day, tot_days, trap_num)) %>% 
  group_by(year, season, genus, trap_status) %>% 
  filter(trap_status == "N") %>% # limit to new individuals
  group_by(year, season) %>% 
  count(genus)

rod_year <- rodents %>%
  dplyr::select(-c(date,trap_day, tot_days, trap_num)) %>% 
  group_by(year, season, genus, trap_status) %>% 
  filter(trap_status == "N", species_code != "Bird") %>% # limit to new individuals, no birds
  group_by(year) %>% 
  count(genus) %>% 
  rename("count" = "n")

#write.csv(rod_year, file = "Data/rodents_by_year.csv")

rod_season <- rodents %>%
  dplyr::select(-c(date,trap_day, tot_days, trap_num)) %>% 
  group_by(year, season, genus, trap_status) %>% 
  filter(trap_status == "N", species_code != "Bird") %>% # limit to new individuals, no birds
  group_by(year, season) %>% 
  count(genus) %>% 
  rename("count" = "n")

#write.csv(rod_season, file = "Data/rodents_by_season.csv")

# Survival vs small mammals
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
tot_surv_yr <- seedlings_obs %>% 
  group_by(precip, cohort) %>% 
  summarise(mean_surv = mean(survival),
            sd_surv = sd(survival),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_yr

rod_grouped <- rod_year %>%
  group_by(year) %>% 
  summarise(rod_count = sum(count))

sm_surv <- cbind(rod_grouped, tot_surv_yr)

# lin reg vs log
# across all precip tx
summary(lm(10*mean_surv~(rod_count), data = sm_surv)) #r2 = 0.59
summary(lm(10*mean_surv~log(rod_count), data = sm_surv))# log improves fit, r2 = 0.69  

# mean survival vs trapped small mammals for each year and ppt?
# sm_surv_fig <- sm_surv %>% 
#   mutate(precip = recode_factor(precip,
#                                 "Control" = "Ambient",
#                                 "IR" = "Wet",
#                                 "RO" = "Drought", .ordered = TRUE)) %>% 
#   ggplot(mapping = aes(x = rod_count, y = 10*mean_surv)) +
#   geom_point(mapping = aes(color = year, shape = precip), size = 10) +
#   geom_smooth(method = "glm", formula = y~log(x), se = F, color = "black", size = 2) +
#   scale_color_manual(values = c("#b0e8f5", "#0033FF", "#169cf0")) + # blue gets darker for the most monsoon rainfall
#   labs(y = "Mean Survival (%)",
#        x = "Trapped",
#        color = "Year",
#        shape = "PPTx") +
#   xlim(0, 250)+
#   ylim(0, 40)+
#   theme_pubr(legend = c("right"))+
#   labs_pubr(base_size = 24) 
# 
# sm_surv_fig

sm_surv_fig_simple <- sm_surv %>% 
  mutate(precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = TRUE)) %>% 
  ggplot(mapping = aes(x = rod_count, y = (10*mean_surv), color = precip))+
  geom_point(size = 6)+
  scale_color_manual(values = c("grey30","blue1","#ba7525")) +
  geom_smooth(method = "glm", formula = y ~ log(x), se = F, size = 2)+
  labs(y = "Mean Survival (%)",
       x = "Trapped",
       color = "PPTx") +
  xlim(0, 250) +
  ylim(0, 40) +
  theme_pubr(legend = c("right"))+
  labs_pubr(base_size = 24)

sm_surv_fig_simple

ggsave(filename = "Figures_Tables/line_sm_surv.tiff",
       plot = sm_surv_fig_simple,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# USDA Rainout - Santa Rita Experimental Range
# Ants Tables
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2021-01-13

# load packages
library(tidyverse)
library(vroom)
library(gt)
library(glue)
library(ggpubr)
library(car)

# read in data
ants <- vroom("Data/ants_combined.csv",
                 col_types = c(.default = "f",
                               samp_date = "D",
                               number = "i"))

glimpse(ants)

# clean data to only be monsoon (comparable across years), remove 2016 and missing data
ants_new <- ants %>%
  dplyr::select(-c(samp_date,sample, order, species)) %>%
  filter(season == "monsoon") %>% 
  group_by(year, season, tx, clip, genus) %>% 
  filter(genus != "",
         year != "2016") %>%
  drop_na(number) %>% 
  summarise(total = sum(number))

# look at summary of cleaned data by ppt treatments
ants_new_tx <- ants_new %>% 
  group_by(year, tx, clip, genus) %>% 
  summarise(total_tx = sum(total))

summary(ants_new)

# look at summary of all ants (no genus) by clip and ppt treatments
ants_total <- ants_new %>%
  group_by(year, clip, tx) %>% 
  summarise(total = sum(total))

# compare clip treatment number of ants
ants_total %>% 
  group_by(clip) %>% 
  summarise(sum(total))

# compare ppt treatment number of ants
ants_total %>% 
  group_by(tx) %>% 
  summarise(sum(total))

# compare number of ants b/w years
ants_total %>% 
  group_by(year) %>% 
  summarise(sum(total))

# differences are not large, are data normal?
hist(ants_new$total) #needs transformation
hist(log(ants_new$total)) #looks more normal

shapiro.test(log(ants_new$total)) # sig. p-value, not normal

kruskal.test(log(total)~year, data = ants_new) #not sig
kruskal.test(log(total)~clip, data = ants_new) #not sig
kruskal.test(log(total)~tx, data = ants_new) #not sig

#clipping not sig so show totals by genus in table by year and ppt treatment for easier viewing

#write.csv(ants_new_tx, file = "Data/ants_by_year_tx.csv")

# Survival vs ants
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

# create data set for survival vs ants fig
tot_surv_yr <- seedlings_obs %>% 
  mutate(year = dplyr::recode(cohort,
                "1" = "2017",
                "2" = "2018",
                "3" = "2019"),
         precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = F)) %>% 
  group_by(precip, clip, year) %>% 
  summarise(mean_surv = mean(survival),
            sd_surv = sd(survival),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_yr

# match col names of ant data to seedling data
ants_total <- ants_total %>% dplyr::rename(precip = tx, total_ants = total)

surv_ants <- full_join(tot_surv_yr, ants_total)

surv_ants

# create fig of just precip and year
hist(surv_ants$total_ants)
hist((surv_ants$mean_surv))

seed_surv <- seedlings_obs %>% 
  mutate(year = dplyr::recode(cohort,
                              "1" = "2017",
                              "2" = "2018",
                              "3" = "2019"),
         precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = F)) %>% 
  group_by(precip, year) %>% 
  summarise(mean_surv = mean(survival),
            sd_surv = sd(survival),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

seed_surv

# calc ants
ants_total_gp <- ants_new %>%
  group_by(year, tx) %>% 
  summarise(total = sum(total))

# rename cols to match seedling data
ants_total_gp <- ants_total_gp %>% dplyr::rename(precip = tx, total_ants = total)

# join seedling and ant data together
surv_ants_gp <- full_join(seed_surv, ants_total_gp)

surv_ants_gp

# create fig for just precip and year
surv_ants_gp_fig <- surv_ants_gp %>%
  ggplot(mapping = aes(x = total_ants, y = (10*mean_surv), color = precip))+
  geom_point(position = "jitter", size = 8,
             aes(shape = year))+
  scale_color_manual(values = c("grey30","blue1","#ba7525")) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Ants Captured",
       color = "PPTx",
       shape = "Year") +
  xlim(0, NA) +
  ylim(0, 40) +
  theme_pubr(legend = c("right"))+
  labs_pubr(base_size = 24)

surv_ants_gp_fig

ggsave(filename = "Figures_Tables/line_ants_surv_grouped.tiff",
       plot = surv_ants_gp_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# with clipping ant data added
surv_ants_fig <- surv_ants %>% 
  ggplot(mapping = aes(x = total_ants, y = (10*mean_surv),color = precip)) +
  geom_point(position = "jitter", size = 8,
             aes(alpha = clip,
                 shape = year))+
  scale_alpha_manual(values = c(0.7, 1))+
  scale_color_manual(values = c("grey30","blue1","#ba7525")) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Ants Captured",
       color = "PPTx",
       shape = "Year",
       alpha = "Clipping") +
  xlim(0, NA) +
  ylim(0, 40) +
  theme_pubr(legend = c("right"))+
  labs_pubr(base_size = 24)

surv_ants_fig

ggsave(filename = "Figures_Tables/line_ants_surv.tiff",
       plot = surv_ants_fig,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

# stats
# across all precip tx
hist(surv_ants$total_ants) # looks normal
summary(lm(mean_surv~(total_ants)*year, data = surv_ants)) #r2 = 0.73
summary(lm(mean_surv~log(total_ants)*year, data = surv_ants))# log improves fit, r2 = 0.75  

# diagnostic plots
ant_mod <- lm(mean_surv~log(total_ants)*year, data = surv_ants)
plot(ant_mod)

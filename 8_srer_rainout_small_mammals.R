# USDA Rainout - Santa Rita Experimental Range
# Small Mammal Trapping Tables
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-12-10

# load packages
library(tidyverse)
library(lubridate)
library(vroom)
library(ggpubr)
library(gt)
library(glue)
library(ggpmisc)
library(agricolae)

# read in data
rodents <- vroom("Data/Trapping_Combined_2017-2019.csv",
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
seedlings <- vroom("Data/seedlings_combined.csv",
                   col_select = -c(1),
                   col_types = c(.default = "f",
                                 date = "D"))
str(seedlings)

# Fate counts
seedlings_fate_full <- seedlings %>% 
  group_by(block, precip, clip, excl, date, cohort) %>%
  count(fate) %>% 
  pivot_wider(names_from = fate,
              values_from = n,
              values_fill = 0) %>% 
  rename(no_germ = "0",
         survival = "1",
         died = "2") %>% 
  mutate(tot_germination = survival + died,
         surv_perc = (survival/tot_germination)) %>% 
  mutate(surv_perc = replace_na(surv_perc, 0))


seedling_fate_post_co1y1 <- seedlings_fate_full %>% 
  filter(cohort == 1) %>% 
  filter(date > "2017-07-21" & date < "2018-07-10") %>%  # emerge (2 weeks following planting) to planting of cohort 2
  mutate(year = as.factor(1))

seedling_fate_post_co1y2 <- seedlings_fate_full %>% 
  filter(cohort == 1) %>% 
  filter(date > "2018-07-10" & date < "2019-08-01") %>%  # planting of cohort 2 to planting of cohort 3
  mutate(year = as.factor(2))

seedling_fate_post_co1y3 <- seedlings_fate_full %>% 
  filter(cohort == 1) %>% 
  filter(date > "2019-08-01") %>% # planting of cohort 3 to end
  mutate(year = as.factor(3))

seedling_fate_post_co2y1 <- seedlings_fate_full %>% 
  filter(cohort == 2) %>% 
  filter(date > "2018-07-26" & date < "2019-08-01") %>%  # emerge to planting of cohort 3
  mutate(year = as.factor(1))

seedling_fate_post_co2y2 <- seedlings_fate_full %>% 
  filter(cohort == 2) %>% 
  filter(date > "2019-08-01") %>%  # planting of cohort three to end
  mutate(year = as.factor(2))

seedling_fate_post_co3y1 <- seedlings_fate_full %>% 
  filter(cohort == 3) %>% 
  filter(date > "2019-08-17") %>%  # emerge of cohort 3 to end
  mutate(year = as.factor(1))

seedling_fate_year <- rbind(seedling_fate_post_co1y1, seedling_fate_post_co1y2,
                            seedling_fate_post_co1y3, seedling_fate_post_co2y1,
                            seedling_fate_post_co2y2, seedling_fate_post_co3y1)


# create data set for raw survival data
tot_surv_yr <- seedling_fate_year %>%
  mutate(cohort = dplyr::recode(cohort,
                              "1" = "2017",
                              "2" = "2018",
                              "3" = "2019"),
         precip = recode_factor(precip,
                                "Control" = "Ambient",
                                "IR" = "Wet",
                                "RO" = "Drought", .ordered = F)) %>% 
  group_by(precip, cohort, year) %>% 
  summarise(mean_surv = mean(surv_perc),
            sd_surv = sd(surv_perc),
            counts = n(),
            se_surv = (sd_surv/sqrt(counts))) %>%
  mutate(upper = mean_surv + se_surv,
         lower = mean_surv - se_surv)

tot_surv_yr

rod_grouped <- rod_year %>%
  group_by(year) %>% 
  summarise(rod_count = sum(count)) %>% 
  rename(cohort = year)

a<-tot_surv_yr %>% filter(year == 1) %>% left_join(rod_grouped)

b<-tot_surv_yr %>% filter(year == 2 & cohort == 2017) %>% add_column(rod_count = c(112))

c<-tot_surv_yr %>% filter(year == 3) %>% add_column(rod_count = c(238))

d<-tot_surv_yr %>% filter(year == 2 & cohort == 2018) %>% add_column(rod_count = c(238))

sm_surv <- rbind(a,b,c,d)

# lin reg vs log
# across all precip tx
glimpse(sm_surv)

hist((sm_surv$rod_count))
hist(log(sm_surv$rod_count)) # doesn't improve but follows relationship between ants/survival


# include cohort or not?
rod_mod <- lm(mean_surv~log(rod_count)+cohort, data =sm_surv) # log improves fit, doesn't change relationships

summary(rod_mod)

rod_mod_noco <- lm(mean_surv~log(rod_count), data =sm_surv) # log improves fit, doesn't change relationships

summary(rod_mod_noco) # include on fig

anova(rod_mod, rod_mod_noco)

# no sig difference to include cohort, going with "noco"/no cohort since simpler model

# insig outside of 1st year survival
rod_mod_yearsurv <- lm(mean_surv~log(rod_count)+year+cohort, data = sm_surv)

summary(rod_mod_yearsurv )

# no sig diff between all years of survival beside 2017
post_rod <- HSD.test(rod_mod_yearsurv, "cohort")
post_rod

# no sig between years planted
post_rod_simple <- HSD.test(rod_mod , "cohort")
post_rod_simple

sm_surv_fig_simple <- sm_surv %>% 
  filter(year == 1) %>% 
  ggplot(mapping = aes(x = rod_count, y = (100*mean_surv), color = precip)) +
  geom_point(size = 8, aes(shape = cohort), position = "jitter")+
  scale_color_manual(values = c("grey30","blue1","#ba7525")) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Rodents Captured",
       color = "PPTx",
       shape = "Year") +
  xlim(0,300) +
  ylim(0, 70) +
  # ggpmisc::stat_poly_eq(formula = y ~ log(x), 
  #                       aes(label =  paste(stat(eq.label),
  #                                          stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
  #                       parse = TRUE)+
  theme_pubr(legend = c("right"))+
  labs_pubr(base_size = 24)

sm_surv_fig_simple

ggsave(filename = "Figures_Tables/line_sm_surv_fig5b.tiff",
       plot = sm_surv_fig_simple,
       dpi = 800,
       width = 22,
       height = 12,
       units = "in",
       compression = "lzw")

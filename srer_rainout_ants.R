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
library(ggpmisc)
library(agricolae)

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
  summarise(total = sum(total)) %>% 
  rename(cohort = year)

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
  group_by(cohort) %>% 
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

# match col names of ant data to seedling data
ants_total <- ants_total %>% dplyr::rename(precip = tx, total_ants = total)

surv_ants <- full_join(tot_surv_yr, ants_total)

surv_ants

# create fig of just precip and year
hist(surv_ants$total_ants)
hist((surv_ants$mean_surv))

surv_ants_gp <- surv_ants %>% 
  group_by(precip, cohort, mean_surv, year) %>% 
  summarise(total_ants = sum(total_ants))

# create fig for just precip and year
surv_ants_gp_fig <- surv_ants_gp %>%
  group_by(precip, cohort) %>% 
  filter(year == 1) %>% 
  #filter(precip != "Drought") %>% 
  ggplot(mapping = aes(x = total_ants, y = (100*mean_surv), color = precip))+
  geom_point(position = "jitter", size = 8,
             aes(shape = cohort))+
  scale_color_manual(values = c("grey30","blue1","#ba7525")) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = F, size = 2)+
  labs(y = "Seedling Survival (%)",
       x = "Ants Captured",
       color = "PPTx",
       shape = "Year") +
  ggpmisc::stat_poly_eq(formula = y ~ (x),
                        aes(label =  paste(stat(eq.label),
                                           stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
                        parse = TRUE)+
  xlim(0, 2500) +
  ylim(0, 70) +
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
summary(lm(mean_surv~(total_ants)+cohort, data = surv_ants)) #r2 = 0.24, fits drought poorly
summary(lm(mean_surv~log(total_ants)+cohort, data = surv_ants))# log improves fit, r2 = 0.24  

# diagnostic plots
ant_mod <- lm(mean_surv~log(total_ants)+cohort, data = surv_ants)
summary(ant_mod)


# insig outside of 1st year survival
ant_mod_2 <- lm(mean_surv~log(total_ants)*year, data = surv_ants)
summary(ant_mod_2)

post_ant <- HSD.test(ant_mod_2, "year")
post_ant

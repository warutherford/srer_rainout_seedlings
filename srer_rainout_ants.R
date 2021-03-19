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
  group_by(year, tx, genus) %>% 
  summarise(total_tx = sum(total))

summary(ants_new)

# look at summary of all ants (no genus) by clip and ppt treatments
ants_total <- ants_new %>%
  group_by(year, tx, clip) %>% 
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

# interactions significant?
summary(glm(log(total)~year*tx*clip, data = ants_new)) #not sig
summary(aov(log(total)~year*tx*clip, data = ants_new)) #not sig

#clipping not sig so show totals by genus in table by year and ppt treatment for easier viewing

#write.csv(ants_new_tx, file = "Data/ants_by_year_tx.csv")

ants_new %>% 
  select(year, total, tx) %>% 
  group_by(year, tx) %>% 
  summarise(total_tx = sum(total),
            # for se, 3 ppt x 3 years
            se = (sd(total)/sqrt(9)))

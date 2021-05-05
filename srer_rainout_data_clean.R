# USDA Rainout - Santa Rita Experimental Range
# Data Cleaning for all three seedling cohorts
# Austin Rutherford
# email:arutherford@email.arizona.edu
# 2020-09-21

# Load packages
library(tidyverse)
library(vroom)

### Read in cohort 1 (seedlings from 2017), make all columns factors
seedlings_1 <- vroom("Data/seedlings_2017.csv",  col_types = c(.default = "f"))

summary(seedlings_1)
str(seedlings_1)

seedlings_1_clean <- seedlings_1 %>%
  # select all columns except date of most recent collection
  select(-c("2020-08-29")) %>% 
  # convert wide data format to long, excluding all but date columns,
  # creating a new "date" column from headers and new "recruit" column of values
  pivot_longer(cols = -c(block, precip, clip, excl, side, rep),
               names_to = "date",
               values_to = "recruit") %>% 
  # Make "date" column a date vector with the ISO8601 format
  mutate(date = lubridate::as_date(date, format = '%Y-%m-%d'),
         # recode recruit values to fix data entry typos and make the dead
         # seedlings to be coded as "2" (recruited = 1, no germination = 0)
         # X in HerbX means seedling had herbivory occur and then died
         recruit = recode(recruit,
                          "1-herbX" = "1-HerbX",
                          "1-Herbx" = "1-HerbX",
                          "Died" = "2")) %>%
  # separate recruitment code into two columns, "fate" and "herbivory" using "-"
  # as the separator
  separate(recruit, c("fate", "herbivory"), "-", fill = "warn") %>% 
  # make herbivory a factor and replace NAs from above with "0" noting no
  # herbivory occured
  mutate(herbivory = as.factor(replace_na(herbivory, 0)),
         # recode herbivory values so 0 = no herbivory, 1 = herbivory,
         # 2 = herbivory and died from it
         herbivory = recode(herbivory,
                            "Herb" = "1",
                            "Herb(lvs)" = "1",
                            "Herb(lvs)(cts)" = "1",
                            "HerbX" = "2",
                            "Herb(lvs)X" = "2",
                            "Herb(lvs)(cts)X" = "2"),
         # create a new column called granivory where "fate" has a value of
         # "Herb" and set to a new value of 1 (NB granivory or consumed seeds
         # was marked only as "Herb" in the field)
         granivory = case_when(fate == 'Herb' ~ '1'),
         # make granivory a factor and replace NAs from above as 0's
         granivory = as.factor(replace_na(granivory, 0)),
         # make fate a factor and recode "Herb" values to 0 since these were
         # noting seeds that were consumed (e.g., granivory)
         fate = as.factor(recode(fate, "Herb" = "0")),
         # create new "cohort" variable with 1 for 1st cohort planted in 2017
         cohort = factor("1"))

str(seedlings_1_clean)

### Read in cohort 2 (seedlings from 2018), make all columns factors
seedlings_2 <- vroom("Data/seedlings_2018.csv",  col_types = c(.default = "f"))

str(seedlings_2)

seedlings_2_clean <- seedlings_2 %>%
  # select all columns except date of most recent collection
  select(-c("2020-08-29")) %>% 
  # convert wide data format to long, excluding all but date columns,
  # creating a new "date" column from headers and new "recruit" column of values
  pivot_longer(cols = -c(block, precip, clip, excl, side, rep),
               names_to = "date",
               values_to = "recruit") %>% 
  # Make "date" column a date vector with the ISO8601 format
  mutate(date = lubridate::as_date(date, format = '%Y-%m-%d'),
         # recode recruit values to fix data entry typos and make the dead
         # seedlings to be coded as "2" (recruited = 1, no germination = 0)
         recruit = recode(recruit,
                          "DIed" = "2",
                          "Died" = "2")) %>%
  # separate recruitment code into two columns, "fate" and "herbivory" using "-"
  # as the separator
  separate(recruit, c("fate", "herbivory"), "-", fill = "warn") %>% 
  # make herbivory a factor and replace NAs from above with "0" noting no
  # herbivory occured
  mutate(herbivory = as.factor(replace_na(herbivory, 0)),
         # recode herbivory values so 0 = no herbivory, 1 = herbivory,
         # 2 = herbivory and died from it 
         herbivory = recode(herbivory,
                            "Herb" = "1",
                            "Herb(lvs)" = "1",
                            "Herb(lvs)(cts)" = "1",
                            "HerbX" = "2",
                            "Herb(lvs)X" = "2",
                            "Herb(lvs)(cts)X" = "2"),
         # create a new column called granivory where "fate" has a value of
         # "Herb" and set to a new value of 1 (NB granivory or consumed seeds
         # was marked only as "Herb" in the field)
         granivory = case_when(fate == 'Herb' ~ '1'),
         # make granivory a factor and replace NAs from above as 0's
         granivory = as.factor(replace_na(granivory, 0)),
         # make fate a factor and recode "Herb" values to 0 since these were
         # noting seeds that were consumed (e.g., granivory)
         fate = as.factor(recode(fate, "Herb" = "0")),
         # create new "cohort" variable with 2 for 2nd cohort planted in 2018
         cohort = factor("2"))

str(seedlings_2_clean)

### Read in cohort 3 (seedlings from 2019), make all columns factors
seedlings_3 <- vroom("Data/seedlings_2019.csv", col_select = c(1:23),
                     col_types = c(.default = "f"))

str(seedlings_3)

seedlings_3_clean <- seedlings_3 %>%
  # select all columns except date of most recent collection
  select(-c("2020-08-29")) %>% 
  # convert wide data format to long, excluding all but date columns,
  # creating a new "date" column from headers and new "recruit" column of values
  pivot_longer(cols = -c(block, precip, clip, excl, side, rep),
               names_to = "date",
               values_to = "recruit") %>% 
  # Make "date" column a date vector with the ISO8601 format
  mutate(date = lubridate::as_date(date, format = '%Y-%m-%d'),
         # recode recruit values to fix data entry typos and make the dead
         # seedlings to be coded as "2" (recruited = 1, no germination = 0)
         recruit = recode(recruit,
                          "Died" = "2")) %>%
  # separate recruitment code into two columns, "fate" and "herbivory" using "-"
  # as the separator
  separate(recruit, c("fate", "herbivory"), "-", fill = "warn") %>% 
  # make herbivory a factor and replace NAs from above with "0" noting no
  # herbivory occured
  mutate(herbivory = as.factor(replace_na(herbivory, 0)),
         # recode herbivory values so 0 = no herbivory, 1 = herbivory,
         # 2 = herbivory and died from it
         herbivory = recode(herbivory,
                            "Herb" = "1",
                            "Herb(lvs)" = "1",
                            "Herb(lvs)(cts)" = "1",
                            "HerbX" = "2",
                            "Herb(lvs)X" = "2",
                            "Herb(lvs)(cts)X" = "2"),
         # create a new column called granivory where "fate" has a value of
         # "Herb" and set to a new value of 1 (NB granivory or consumed seeds
         # was marked only as "Herb" in the field)
         granivory = case_when(fate == 'Herb' ~ '1'),
         # make granivory a factor and replace NAs from above as 0's
         granivory = as.factor(replace_na(granivory, 0)),
         # make fate a factor and recode "Herb" values to 0 since these were
         # noting seeds that were consumed (e.g., granivory)
         fate = as.factor(recode(fate, "Herb" = "0")),
         # create new "cohort" variable with 3 for 3rd cohort planted in 2019
         cohort = factor("3"))

str(seedlings_3_clean)

# Combine three seedling cohort data sets together (append rows)
seedlings <- bind_rows(seedlings_1_clean, seedlings_2_clean, seedlings_3_clean)

# Save combined data set as new csv
#write.csv(seedlings, file = "Data/seedlings_combined.csv")

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

# not grouping...why?
rod_tbl_year <- rod_year %>% 
  drop_na(genus) %>% 
  gt(rowname_col = "Genera",
     groupname_col = "year") %>%
  summary_rows(groups = "year",
               columns = vars(count),
               fns = list(total = "sum")) %>% 
  tab_header(title = "Small Mammal Surveys") %>% 
  cols_label(count = "New Captures",
             genus = "Genera") 

rod_tbl_season <- rod_year %>% 
  gt(rowname_col = "Genera",
     groupname_col = c("year", "season")) %>%
  tab_header(title = "Small Mammal Surveys") %>% 
  cols_label(count = "New Captures",
             genus = "Genera")

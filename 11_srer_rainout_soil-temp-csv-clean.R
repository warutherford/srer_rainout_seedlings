# The credit for the script below belongs to Bryan Urban https://gist.github.com/bjurban/ef74a5accf42c43480a1
# Combine multiple raw HOBO data files into a single .csv file 
# Author: Bryan Urban
# Email: burban@fraunhofer.org
# Date: 2015-04-30

## SETUP --------------------

library(data.table)
library(tidyverse)
library(lubridate)
library(udunits2)

# change these to match the folders containing the data
raw_dir <- "./Data/site-env-data/soil-temp/raw"
out_dir <- "./Data/site-env-data/soil-temp" 


## LOAD RAW DATA ------------
# get file names: 
pattern = ".*csv$" # for identifying files to read
fns <- list.files(raw_dir, pattern=pattern, full.names=TRUE)

# load data into lists
read_and_label <- function(x,...){
  z <- fread(x,...)
  
  # add file name without the extension as id column
  pattern <- "(.*\\/)([^.]+)(\\.csv$)"
  z$ids <- sub(pattern, "\\2", x)
  z
}

# reads columns 2 and 3 (timestamp and temperature) into a list of data.table
all_data <- 
  lapply(fns, function(x,...) {read_and_label(x,...)},
          header = TRUE, skip=0
  )

## PROCESS RAW DATA ---------
# drop errors, merge into one large data.table, make table a data frame, name columns
all_data <- all_data[sapply(all_data, is.data.frame)]
all_data <- rbindlist(all_data, fill = TRUE, use.names = TRUE, idcol = TRUE)
all_data <- as.data.frame(all_data)

colnames(all_data)
setnames(all_data, c("setID","datetime","soiltempC", "ids", "lux"))

# transform to long format, force date format, make column for treatments, make treatments factors
comb_data_amb <- all_data %>% 
  mutate(datetime = ymd_hms(datetime)) %>% 
  separate("ids", c("precip","location", "side"), sep = "-", remove = TRUE, fill = "warn") %>%
  separate("side", c("section", "extra", "extraX"), sep = "_", fill = "right") %>% # 'extra' and 'extraX' are dummy columns from the file names
  mutate(precip = as.factor(precip),
         clip = as.factor(section),
         location = as.factor(location)) %>% 
  mutate(precip = recode_factor(precip,
                         "CO" = "Ambient",
                         "RO" = "Drought",
                         "IR" = "Wet")) %>% 
  filter(precip == "Ambient", .preserve = TRUE) %>% 
  mutate(clip = recode_factor(section, 
                       "E1" = "Clipped",
                       "E2" = "Clipped",
                       "W1" = "Unclipped",
                       "W2" = "Unclipped")) %>% 
  dplyr::select(-c(extra, extraX)) # remove dummy columns from file names

comb_data_d_w <- all_data %>% 
  mutate(datetime = ymd_hms(datetime)) %>% 
  separate("ids", c("precip","location", "side"), sep = "-", remove = TRUE, fill = "warn") %>%
  separate("side", c("section", "extra", "extraX"), sep = "_", fill = "right") %>% # 'extra' and 'extraX' are dummy columns from the file names
  mutate(precip = as.factor(precip),
         clip = as.factor(section),
         location = as.factor(location)) %>% 
  mutate(precip = recode_factor(precip,
                                "CO" = "Ambient",
                                "RO" = "Drought",
                                "IR" = "Wet")) %>% 
  filter(precip == "Drought" | precip == "Wet") %>% 
  mutate(clip = recode_factor(section, 
                              "E1" = "Unclipped",
                              "E2" = "Unclipped",
                              "W1" = "Clipped",
                              "W2" = "Clipped")) %>% 
  dplyr::select(-c(extra, extraX)) # remove dummy columns from file names
  
comb_data  <- full_join(comb_data_amb, comb_data_d_w)

## SAVE PROCESSED DATA ------
# write data in one big file
write.csv(comb_data, paste(out_dir, "all_soil_temp_light.csv", sep="/"), 
          row.names=FALSE)

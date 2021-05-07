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
setnames(all_data, c("setID","datetime","soiltempC", "ids"))

# transform to long format, force date format, make column for treatments, make treatments factors
comb_data <- all_data %>% 
  mutate(datetime = ymd_hms(datetime)) %>% 
  separate("ids", c("precip","location", "side"), sep = "-", remove = TRUE, fill = "warn") %>%
  separate("side", c("section", "extra"), sep = "_", fill = "right") %>% 
  mutate(precip = as.factor(precip),
         clip = as.factor(section),
         location = as.factor(location)) %>% 
  mutate(precip = recode_factor(precip,
                         "CO" = "Ambient"),
                         #"RO" = "Drought",
                         #"IR" = "Wet"),
         clip = recode_factor(section, 
                       "E1" = "Clipped",
                       "E2" = "Clipped",
                       "W1" = "Unclipped",
                       "W2" = "Unclipped")) %>% 
  dplyr::select(-extra)
  

## SAVE PROCESSED DATA ------
# write data in one big file
write.csv(comb_data, paste(out_dir, "all_soil_temp.csv", sep="/"), 
          row.names=FALSE)

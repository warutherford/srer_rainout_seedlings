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
raw_dir <- "./Data/site-env-data/soil-moisture/raw"
out_dir <- "./Data/site-env-data/soil-moisture" 


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
  lapply(fns, function(x,...) {try(read_and_label(x,...))},
         select=1:5, header=TRUE, skip=0
  ) 

## PROCESS RAW DATA ---------
# drop errors, merge into one large data.table, make table a data frame, name columns
all_data <- all_data[sapply(all_data, is.data.frame)]
all_data <- rbindlist(all_data, fill = TRUE, use.names = TRUE, idcol = TRUE)
all_data <- as.data.frame(all_data)

colnames(all_data)
setnames(all_data, c("setID","datetime","Control_Unclipped",
                     "Control_Clipped", "Control_Clipped", "Control_Unclipped",
                     "ids", "Control_Unclipped", "Control_Clipped", "Control_Clipped",
                     "Control_Clipped", "IR_Clipped", "IR_Clipped", "IR_Unclipped", "IR_Unclipped",
                     "IR_Clipped", "RO_Clipped", "RO_Clipped", "RO_Unclipped", "RO_Unclipped", "RO_Unclipped"))

# transform to long format, force date format, make column for treatments, make treatments factors
comb_data <- all_data %>% 
  pivot_longer(names_to = "treatment",
               values_to = "moisture",
               values_drop_na = TRUE,
               -c(setID, datetime, ids)) %>% 
  mutate(datetime = parse_date_time(datetime, "%y-%m-%d %H:%M:%S")) %>% 
  separate("treatment", c("precip", "side"), sep = "_", remove = TRUE, fill = "warn") %>% 
  mutate(precip = as.factor(precip),
         clip = as.factor(side)) %>% 
  dplyr::select(datetime, precip, clip, moisture, ids)

# split into two
comb_data_1 <- comb_data %>% 
  slice(1:1124168)

comb_data_2 <- comb_data %>% 
  slice(1124169:2248336)

## SAVE PROCESSED DATA ------
# write data in one big file for use later
# write.csv(comb_data, paste(out_dir, "all_sm.csv", sep="/"), 
#           row.names=FALSE)

# write data in 2 files for upload
write.csv(comb_data_1, paste(out_dir, "all_sm_1.csv", sep="/"), 
          row.names=FALSE)

write.csv(comb_data_2, paste(out_dir, "all_sm_2.csv", sep="/"), 
          row.names=FALSE)

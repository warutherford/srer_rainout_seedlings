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
raw_dir <- "./Data/site-env-data/ppt-temp/raw"
out_dir <- "./Data/site-env-data/ppt-temp" 


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
         select=1:3, header=FALSE, skip=1
  ) 

## PROCESS RAW DATA ---------
# drop errors, merge into one large data.table, make table a data frame, name columns, parse timestamp
all_data <- all_data[sapply(all_data, is.data.frame)]
all_data <- rbindlist(all_data)
all_data <- as.data.frame(all_data)

setnames(all_data, c("datetime", "tempF","pptmm", "ids"))

comb_data <- all_data %>% 
  mutate(datetime = ymd_hms(datetime),
         tempC = ud.convert(tempF, "degree_fahrenheit", "celsius")) %>% 
  dplyr::select(datetime, tempC, pptmm, ids)


## SAVE PROCESSED DATA ------
# write data in one big file
write.csv(comb_data, paste(out_dir, "all_air_ppt_temp.csv", sep="/"), 
          row.names=FALSE)

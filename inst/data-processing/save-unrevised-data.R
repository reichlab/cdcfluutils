library(dplyr)
library(tidyr)
library(purrr)
library(MMWRweek)
library(lubridate)
source("https://raw.githubusercontent.com/cmu-delphi/delphi-epidata/master/src/client/delphi_epidata.R")

# Fetch ILI data for national level and all HHS regions
location_codes <- c("nat", paste0("hhs", 1:10))

# all issues for year 1997 to present
current_year <- lubridate::year(Sys.time())
all_issues <- expand.grid(
    year = 1997:current_year,
    week = sprintf("%02d", 1:53)
  ) %>%
  apply(1, function(x) paste(x, collapse = "")) %>%
  as.integer() %>%
  sort()

# all epidemic weeks from 199740 to present
temp <- MMWRweek::MMWRweek(Sys.time())
current_epiweek <- paste0(temp$MMWRyear, temp$MMWRweek) %>% as.numeric()
epiweeks_range <- c(199740, current_epiweek)

# fetch the data
flu_data_with_backfill <- cdcfluutils::fetch_delphi_data_multi_issue(
  source = "fluview",
  regions = location_codes,
  issues = all_issues,
  epiweeks_range = epiweeks_range)

# save in package data folder
save(flu_data_with_backfill, file = "data/flu_data_with_backfill.rdata")

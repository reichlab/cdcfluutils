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
nat_reg_flu_data_with_backfill <- cdcfluutils::fetch_delphi_data_multi_issue(
  source = "fluview",
  regions = location_codes,
  issues = all_issues,
  epiweeks_range = epiweeks_range)

# save in package data folder
save(nat_reg_flu_data_with_backfill, file = "data/nat_reg_flu_data_with_backfill.rdata")


# Fetch ILI data for state and local level
location_codes <- c(
  'al', 'ak', 'az', 'ar', 'ca', 'co', 'ct', 'de', 'fl', 'ga', 'hi', 'id', 'il',
  'in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt',
  'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or',
  'pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy',
  'as', 'mp', 'dc', 'gu', 'pr', 'vi', 'ord', 'lax', 'jfk')

# fetch the data
state_local_flu_data_with_backfill <- cdcfluutils::fetch_delphi_data_multi_issue(
  source = "fluview",
  regions = location_codes,
  issues = all_issues,
  epiweeks_range = epiweeks_range)

# save in package data folder
save(state_local_flu_data_with_backfill, file = "data/state_local_flu_data_with_backfill.rdata")



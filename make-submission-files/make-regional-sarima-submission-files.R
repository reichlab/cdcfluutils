library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(forecast)
library(MMWRweek)
library(cdcfluutils)
#library(cdcflu20192020)
library(predx)
#library(FluSight)

### Get command line arguments
#args <- c("FALSE", "none", "1")
args <- commandArgs(trailingOnly = TRUE)
seasonal_difference <- as.logical(args[1])
backfill_method <- args[2]
replicate_ind <- as.integer(args[3])

#analysis_time_year <- "2010"
#analysis_time_week <- "40"
#seasonal_difference <- FALSE
#backfill_method <- "none"
#backfill_method <- "forecast_input"
#backfill_method <- "post-hoc"
#replicate_ind <- 1

method <- paste0("sarima_seasonal_difference_", seasonal_difference,
  "_backfill_", gsub("-", "_", backfill_method))

submissions_save_paths <- paste0("inst/submissions/region-", method)

if(backfill_method == "none") {
  submissions_save_paths <- c(submissions_save_paths,
    paste0("~/Documents/research/epi/flu/cdc-flusight-ensemble/model-forecasts/real-time-component-models/ReichLab_sarima_seasonal_difference_",
    seasonal_difference)
  )
}

data <- download_and_preprocess_flu_data()

data <- data %>% mutate(
  epiweek = as.integer(paste0(year, sprintf("%02d", week)))
)

analysis_time_season <- max(data$season)
analysis_time_epiweek <- max(data$epiweek)
analysis_time_year <- substr(as.character(analysis_time_epiweek), 1, 4)
analysis_time_week <- substr(as.character(analysis_time_epiweek), 5, 6)

## Parameters used in simulating trajectories via sarima
simulate_trajectories_sarima_params <- list(
  fits_filepath = file.path("inst",
    "estimation",
    "region-sarima",
    ifelse(seasonal_difference,
      "fits-seasonal-differencing",
      "fits-no-seasonal-differencing")),
  prediction_target_var = "weighted_ili",
  seasonal_difference = seasonal_difference,
  transformation = "box-cox",
  first_test_season = analysis_time_season,
  age = NA,
  regional_switch = "NatRegState"
)

weeks_in_first_season_year <-
  get_num_MMWR_weeks_in_first_season_year(analysis_time_season)



all_regions <- c(paste0("Region ", 1:10), "National")
res_predx <- purrr::map_dfr(all_regions, function(region) {
  cat(region)
  cat(" ")
  cat(Sys.time())
  cat("\n")
  get_submission_one_region_via_trajectory_simulation(
    data = data,
    analysis_time_season = analysis_time_season,
    first_analysis_time_season_week = 10, # == week 40 of year
    last_analysis_time_season_week = weeks_in_first_season_year - 11, # analysis for 33-week season
    region = region,
    prediction_target_var = "weighted_ili",
    incidence_bins = data.frame(
      lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
      upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
    incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
    n_trajectory_sims = 50000,
    simulate_trajectories_function = sample_predictive_trajectories_arima_wrapper,
    simulate_trajectories_params = simulate_trajectories_sarima_params,
    backfill_adjust = backfill_method) #, "forecast input", "post-hoc"
}) %>%
  mutate(
    location = ifelse(
      location == "National",
      "US National",
      paste0("HHS ", location)
    )
  )

options(digits = 22, scipen = 9999)
csv_df <- cdcfluutils::predx_to_submission_df(
  res_predx,
  ew = analysis_time_week,
  year = analysis_time_year,
  team = "Kernel of Truth")

submissions_save_path <- submissions_save_paths[1]
predx_res_file <- file.path(
  submissions_save_path,
  "predx",
  paste0(
    "EW", analysis_time_week,
    "-", analysis_time_year,
    "-ReichLab_", method,
    "-replicate_", replicate_ind,
    ".rds"))
saveRDS(res_predx, predx_res_file)

csv_res_file <- file.path(
  submissions_save_path,
  "csv",
  paste0(
    "EW", analysis_time_week,
    "-", analysis_time_year,
    "-ReichLab_", method,
    ".csv"))
write.csv(
  csv_df,
  csv_res_file,
  row.names = FALSE)

if(length(submissions_save_paths) > 1) {
  submissions_save_path <- submissions_save_paths[2]
  csv_res_file <- file.path(
    submissions_save_path,
    paste0(
      "EW", analysis_time_week,
      "-", analysis_time_year,
      "-ReichLab_", method,
      ".csv"))
  write.csv(
    csv_df,
    csv_res_file,
    row.names = FALSE)
}

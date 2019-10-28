library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(forecast)
library(kcde)
library(copula)
library(MMWRweek)
library(cdcfluutils)
#library(cdcflu20192020)
library(predx)
#library(FluSight)

### Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
backfill_method <- args[1]

#analysis_time_year <- "2018"
#analysis_time_week <- "47"
#seasonal_difference <- FALSE
#backfill_method <- "none"
#backfill_method <- "forecast_input"
#backfill_method <- "post-hoc"

method <- paste0("kcde_backfill_", gsub("-", "_", backfill_method))

submissions_save_paths <- c(
    paste0("inst/submissions/region-", method),
    paste0("~/Documents/research/epi/flu/cdc-flusight-ensemble/model-forecasts/real-time-component-models/ReichLab_",
      method)
  )

data <- download_and_preprocess_flu_data()

data <- data %>% mutate(
  epiweek = as.integer(paste0(year, sprintf("%02d", week)))
)

analysis_time_season <- max(data$season)
analysis_time_epiweek <- max(data$epiweek)
analysis_time_year <- substr(as.character(analysis_time_epiweek), 1, 4)
analysis_time_week <- substr(as.character(analysis_time_epiweek), 5, 6)

## Parameters used in simulating trajectories via kcde
simulate_trajectories_kcde_params <- list(
  n_kcde_sims = 10^5,
  copula_save_path = "inst/estimation/region-kcde/copula-fits",
  estimation_results_path = "inst/estimation/region-kcde/fits",
  max_lag = "1",
  seasonality = TRUE,
  bw_parameterization = "diagonal",
  last_analysis_time_season_week = 41,
  first_test_season = analysis_time_season
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
    simulate_trajectories_function = simulate_trajectories_kcde,
    simulate_trajectories_params = simulate_trajectories_kcde_params,
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

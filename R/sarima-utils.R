## utility functions for SARIMA fits

#' A wrapper around simulate.sarimaTD suitable for use as the
#' \code{simulate_trajectories_function} argument to
#' \code{get_log_scores_via_trajectory_simulation}.
#'
#' @param n_sims number of trajectories to simulate
#' @param max_prediction_horizon how many steps ahead to simulate
#' @param data data set
#' @param region region
#' @param analysis_time_season season in which we're predicting
#' @param analysis_time_season_week week of the season in which we're making our
#'   predictions, using all data up to the analysis time to make predictions for
#'   later time points
#' @param params other parameters.  A list with the following entries:
#'   * fits_filepath = path to a directory where SARIMA model fits are located
#'   * prediction_target_var = string naming variable in data we are predicting
#'   * seasonal_difference = logical specifying whether a seasonal difference
#'       should be computed manually before passing to auto.arima
#'   * transformation = string, either "log", "box-cox", or "none", indicating
#'       type of transformation to do
#'   * first_test_season = string, in format of "2011/2012", specifying first
#'       test season.
#'
#' @return an n_sims by h matrix with simulated values
#'
#' @export
sample_predictive_trajectories_arima_wrapper <- function(
  n_sims,
  max_prediction_horizon,
  data,
  region,
  analysis_time_season,
  analysis_time_season_week,
  params,
  age,
  regional_switch
) {
  
  require(sarimaTD)
  
  ## load SARIMA fit
  reg_str <- gsub(" ", "_", region)
  
  if (regional_switch != "Hosp"){
  fit_filepath <- file.path(
    params$fits_filepath,
    paste0(
      "sarima-",
      reg_str,
      "-seasonal_difference_", params$seasonal_difference,
      "-transformation_", params$transformation,
      "-first_test_season_",
      gsub("/", "_", params$first_test_season),
      ".rds"))
  
  }
  else{
  fit_filepath <- file.path(
    params$fits_filepath,
    paste0(
      "sarima-",
      reg_str,
      "-age-",age,
      "-seasonal_difference_", params$seasonal_difference,
      "-transformation_", params$transformation,
      "-first_test_season_",
      gsub("/", "_", params$first_test_season),
      ".rds"))
  }
  ## If no SARIMA fit, exit early by returning a matrix of NAs
  if(!file.exists(fit_filepath)) {
    warning("no file found for existing fit.")
    return(matrix(NA, nrow = n_sims, ncol = max_prediction_horizon))
  }
  
  sarima_fit <- readRDS(file = fit_filepath)
  
  inc_trajectory_samples <- sarimaTD:::simulate.sarimaTD(
    sarima_fit,
    nsim = n_sims,
    seed = NULL,
    newdata = data[, params$prediction_target_var],
    h = max_prediction_horizon
  )
  
  return(inc_trajectory_samples)
}


#' Estimate SARIMA model using data up to but not including first_test_season
#'
#' @param data regional dataset with structure like regionflu-cleaned
#' @param reg_num region number for estimation
#' @param first_test_season string indicating first test season
#' @param d order of first differencing
#' @param D order of seasonal differencing
#' @param seasonal_difference boolean; take a seasonal difference before passing
#'   to auto.arima?
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", or "none" (default is "none")
#' @param prediction_target_var character specifying column name of modeled 
#'   variable, defaults to "weighted_ili" 
#' @param path path in which to save files
#'
#' @return NULL just saves files
#'
#' @export
fit_region_sarima <- function(
  data,
  region,
  first_test_season,
  d = NA,
  D = NA,
  seasonal_difference = TRUE,
  transformation = c("none", "box-cox", "log"),
  prediction_target_var = "weighted_ili",
  path) {
    
  transformation <- match.arg(transformation)
    
  require(sarimaTD)
  
  ## subset data to be only the region of interest
  data <- data[data$region == region,]

  ## Subset data to do estimation using only data up to (and not including)
  ## first_test_season.  remainder are held out
  first_ind_test_season <- min(which(data$season == first_test_season))
  data <- data[seq_len(first_ind_test_season - 1), , drop = FALSE]

  sarima_fit <- fit_sarima(
    y = data[, prediction_target_var],
    ts_frequency = 52,
    transformation = transformation,
    seasonal_difference = seasonal_difference,
    d = d,
    D = D)

  filename <- paste0(
    path,
    "sarima-",
    gsub(" ", "_", region),
    "-seasonal_difference_", seasonal_difference,
    "-transformation_", transformation,
    "-first_test_season_", gsub("/", "_", first_test_season),
    ".rds")
  saveRDS(sarima_fit, file = filename)
}

fit_hosp_sarima <- function(
  data,
  age,
  first_test_season,
  d = NA,
  D = NA,
  seasonal_difference = TRUE,
  transformation = c("none", "box-cox", "log"),
  prediction_target_var = "weighted_ili",
  path) {
  
  transformation <- match.arg(transformation)
  
  require(sarimaTD)
  
  ## subset data to be only the region of interest
  
  ## Subset data to do estimation using only data up to (and not including)
  ## first_test_season.  remainder are held out
  first_ind_test_season <- min(which(data$season == first_test_season))
  #data <- data[seq_len(first_ind_test_season - 1), , drop = FALSE]
  sarima_fit <- fit_sarima(
    y = data[data$age_label==age,]$weeklyrate,
    ts_frequency = 52,
    transformation = transformation,
    seasonal_difference = seasonal_difference,
    d = d,
    D = D)
  
  filename <- paste0(
    path,
    "sarima-Entire_Network",
    "-age-",age,
    "-seasonal_difference_", seasonal_difference,
    "-transformation_", transformation,
    "-first_test_season_", gsub("/", "_", first_test_season),
    ".rds")
  saveRDS(sarima_fit, file = filename)
}

fit_hosp_sarima <- function(
  data,
  age,
  first_test_season,
  d = NA,
  D = NA,
  seasonal_difference = TRUE,
  transformation = c("none", "box-cox", "log"),
  prediction_target_var = "weighted_ili",
  path) {
  
  transformation <- match.arg(transformation)
  
  require(sarimaTD)
  
  ## subset data to be only the region of interest
  
  ## Subset data to do estimation using only data up to (and not including)
  ## first_test_season.  remainder are held out
  first_ind_test_season <- min(which(data$season == first_test_season))
  #data <- data[seq_len(first_ind_test_season - 1), , drop = FALSE]
  sarima_fit <- fit_sarima(
    y = data[data$age_label==age,]$weeklyrate,
    ts_frequency = 52,
    transformation = transformation,
    seasonal_difference = seasonal_difference,
    d = d,
    D = D)
  
  filename <- paste0(
    path,
    "sarima-Entire_Network",
    "-age-",substr(age,1,5),
    "-seasonal_difference_", seasonal_difference,
    "-transformation_", transformation,
    "-first_test_season_", gsub("/", "_", first_test_season),
    ".rds")
  saveRDS(sarima_fit, file = filename)
}
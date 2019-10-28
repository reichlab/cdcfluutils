#' Utility function to compute onset week based on a trajectory of incidence values
#'
#' @param incidence_trajectory a numeric vector with incidence for each time
#'   point in a season
#' @param baseline the threshold that incidence must cross to count as an onset
#' @param onset_length number of consecutive time points that incidence must
#'   exceed the baseline threshold in order to count as the season onset
#' @param first_season_week number of weeks in year corresponding to the first
#'   week in the season.  For example, our code takes this value to be 31:
#'   a new influenza season starts on the 31st week of each year.
#' @param weeks_in_first_season_year How many MMWR weeks are in the first year
#'   of the season?  For example, in the 2000/2001 season, the first year is
#'   2000.  There were 52 MMWR weeks in 2000.
#'
#' @return the smallest index i such that every entry of
#'   incidence_trajectory[seq(from = i, length = onset_length)]
#'   is >= baseline, if such an index exists.
#'   Otherwise, the character vector "none"
#'
#' @export
get_onset_week <- function(incidence_trajectory,
                           baseline,
                           onset_length,
                           first_season_week = 31,
                           weeks_in_first_season_year) {
  
  exceeded_threshold <- sapply(
    seq_len(length(incidence_trajectory) - onset_length),
    function(start_ind) {
      above_baseline <- incidence_trajectory[seq(from = start_ind, length = onset_length)] >= baseline
      length(above_baseline)>0 &&
        all(above_baseline) &&
        !all(is.na(incidence_trajectory))
    }
  )
  
  if(any(exceeded_threshold, na.rm = TRUE)) {
    season_week <- min(which(exceeded_threshold))
    
    return(season_week)
  } else {
    return("none")
  }
}



#' Compute season onset, peak week, and peak incidence
#'
#' @param data a data frame containing at minimum columns named season,
#'   season_week and a column with some sort of incidence measure
#' @param season the season to look at
#' @param first_CDC_season_week the first week of the season to use for
#'   calculating onset and peak
#' @param last_CDC_season_week the last week of the season to use for
#'   calculating onset and peak
#' @param onset_baseline numeric baseline value for determining season onset
#' @param incidence_var a character string naming the variable in the data
#'   argument containing a measure of incidence, or an integer index
#' @param incidence_bins a data frame with variables lower and upper defining
#'   lower and upper endpoints to use in binning incidence
#' @param incidence_bin_names a character vector with a name for each incidence
#'   bin
#'
#' @return a list with four entries:
#'   1) observed_onset_week, either an integer between first_CDC_season_week
#'     and last_CDC_season_week (inclusive), or "none"
#'   2) observed_peak_week, an integer between first_CDC_season_week and
#'     last_CDC_season_week (inclusive)
#'   3) observed_peak_inc, a numeric with the maximum value of the specified
#'     incidence measure between first_CDC_season_week and last_CDC_season_week
#'   4) observed_peak_inc_bin, character name of incidence bin for peak incidence
#'
#' @export
get_observed_seasonal_quantities <- function(
  data,
  season,
  first_CDC_season_week = 10,
  last_CDC_season_week = 42,
  onset_baseline,
  incidence_var,
  incidence_bins,
  incidence_bin_names
) {
  first_season_ind <- min(which(data$season == season))
  last_season_ind <- max(which(data$season == season))
  
  obs_inc_in_season_leading_trailing_nas <-
    data[seq(from = first_season_ind, to = last_season_ind),
         incidence_var]
  
  ## pad so that we start at season week 1
  if(data$season_week[first_season_ind] != 1) {
    obs_inc_in_season_leading_trailing_nas <- c(
      rep(NA, data$season_week[first_season_ind] - 1),
      obs_inc_in_season_leading_trailing_nas)
  }
  
  ## set values before first analysis time season week or after last
  ## analysis time season week to NA
  ## these are outside of the bounds of the season the CDC wants to look at
  obs_inc_in_season_leading_trailing_nas[
    seq_len(first_CDC_season_week - 1)] <- NA
  if(length(obs_inc_in_season_leading_trailing_nas) >
     last_CDC_season_week) {
    obs_inc_in_season_leading_trailing_nas[
      seq(from = last_CDC_season_week + 1,
          to = length(obs_inc_in_season_leading_trailing_nas))] <- NA
  }
  
  observed_peak_inc <- max(
    obs_inc_in_season_leading_trailing_nas,
    na.rm = TRUE)
  observed_peak_inc_bin <- get_inc_bin(observed_peak_inc, return_character = TRUE)
  
  ## peak week timing is based on rounded values
  rounded_observed_peak_inc <- round(observed_peak_inc, 1)
  rounded_obs_inc_in_season <- sapply(obs_inc_in_season_leading_trailing_nas,
                                      round, digits = 1)
  
  observed_peak_week <-
    which(rounded_obs_inc_in_season == as.numeric(rounded_observed_peak_inc))
  
  weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(season)
  observed_onset_week <- get_onset_week(
    incidence_trajectory = rounded_obs_inc_in_season,
    #    incidence_trajectory = obs_inc_in_season_leading_trailing_nas, # used in stable method
    baseline = onset_baseline,
    onset_length = 3L,
    first_season_week = 31,
    weeks_in_first_season_year = weeks_in_first_season_year
  )
  
  return(list(observed_onset_week = observed_onset_week,
              observed_peak_week = observed_peak_week,
              observed_peak_inc = observed_peak_inc,
              observed_peak_inc_bin = observed_peak_inc_bin
  ))
}



#' Compute season onset, peak week, and peak incidence
#'
#' @param data a data frame containing at minimum columns named season,
#'   season_week and a column with some sort of incidence measure
#' @param season the season to look at
#' @param incidence_var a character string naming the variable in the data
#'   argument containing a measure of incidence, or an integer index
#'
#' @return a list with four entries:
#'   1) observed_onset_week, either an integer between first_CDC_season_week
#'     and last_CDC_season_week (inclusive), or "none"
#'   2) observed_peak_week, an integer between first_CDC_season_week and
#'     last_CDC_season_week (inclusive)
#'   3) observed_peak_inc, a numeric with the maximum value of the specified
#'     incidence measure between first_CDC_season_week and last_CDC_season_week
#'   4) observed_peak_inc_bin, character name of incidence bin for peak incidence
#'
#' @export
get_official_observed_seasonal_quantities <- function(
  data,
  season,
  incidence_var
) {
  require(FluSight)
  
  first_season_ind <- min(which(data$season == season))
  last_season_ind <- max(which(data$season == season))
  
  results <- FluSight::create_truth(fluview = FALSE,
                                    year = substr(season, 1, 4),
                                    weekILI = data,
                                    challenge = "ilinet")
  
  return(list(observed_onset_week = observed_onset_week,
              observed_peak_week = observed_peak_week,
              observed_peak_inc = observed_peak_inc,
              observed_peak_inc_bin = observed_peak_inc_bin
  ))
}



#' Calculate "log scores" for the purpose of the competition -- log[sum_i(p_i)] where p_i is the model's
#' probability of bin i and i runs over some bins adjacent to the bin where the observed quantity was.
#'
#' @param bin_log_probs named numeric vector with log probability of each bin;
#'   names identify the bins
#' @param observed_bin character vector with name(s) of the observed bin(s)
#'   (Note that peak can occur in multiple bins)
#' @param prediction_target
#'
#' @return log score for the given observation
#'
#' @export
compute_competition_log_score <- function(bin_log_probs,
                                          observed_bin,
                                          prediction_target) {
  ## validate probabilities sum to 1 (log of sum = 0) and if not, force them to, with warning
  log_sum_bin_probs <- logspace_sum(bin_log_probs)
  if(log_sum_bin_probs > 2 * .Machine$double.eps) {
    warning(paste(prediction_target, "probabilities do not sum to 1; sum is", exp(log_sum_bin_probs), "automatically adjusting."))
    bin_log_probs <- bin_log_probs - log_sum_bin_probs
  }
  
  ## validate bin names match expected for prediction_target and
  ## observed_bin has appropriate length
  if(identical(prediction_target, "onset_week")) {
    expected_bin_names_52 <- c(as.character(10:42), "none")
    expected_bin_names_53 <- c(as.character(10:43), "none")
    
    if(!identical(length(observed_bin), 1L) || !identical(typeof(observed_bin), "character")) {
      stop("For prediction target onset_week, observed_bin must be a character vector of length 1")
    }
    
    if(identical(sort(expected_bin_names_52), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_52
    } else if(identical(sort(expected_bin_names_53), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_53
    } else {
      stop("invalid names for the vector bin_log_probs")
    }
  } else if(identical(prediction_target, "peak_week")) {
    expected_bin_names_52 <- as.character(10:42)
    expected_bin_names_53 <- as.character(10:43)
    
    if(length(observed_bin) == 0 || !identical(typeof(observed_bin), "character")) {
      stop("For prediction target onset_week, observed_bin must be a character vector of length > 0")
    }
    
    if(identical(sort(expected_bin_names_52), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_52
    } else if(identical(sort(expected_bin_names_53), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_53
    } else {
      stop("invalid names for the vector bin_log_probs")
    }
  } else if(prediction_target %in% c("peak_inc", paste0("ph", 1:4, "_inc"))) {
    expected_bin_names <- as.character(seq(from = 0, to = 13, by = 0.1))
    
    if(!identical(length(observed_bin), 1L) || !identical(typeof(observed_bin), "character")) {
      stop("For given prediction target, observed_bin must be a character vector of length 1")
    }
    
    if(!identical(sort(expected_bin_names), sort(names(bin_log_probs)))) {
      stop("invalid names for the vector bin_log_probs")
    }
  } else {
    stop("Invalid value for prediction_target: must be 'onset_week', 'peak_week', 'peak_inc', or 'phk_inc' for a horizon k")
  }
  
  ## validate observed_bin is one of the expected_bin_names
  if(!(all(observed_bin %in% expected_bin_names))) {
    stop(paste0(
      "observed_bin must be one of (",
      paste(expected_bin_names, collapse = ", "),
      ")"
    ))
  }
  
  ## get bins to sum over
  obs_bin_inds <- sapply(observed_bin, function(bin_name) {
    which(expected_bin_names == bin_name)
  })
  if(identical(prediction_target, "onset_week")) {
    if(identical(observed_bin, "none")) {
      bins_to_sum <- obs_bin_inds
    } else {
      bins_to_sum <- obs_bin_inds + seq(from = -1, to = 1, by = 1)
      bins_to_sum <- bins_to_sum[
        bins_to_sum >= 1 & bins_to_sum <= length(expected_bin_names) - 1]
    }
  } else if(identical(prediction_target, "peak_week")) {
    bins_to_sum <- unique(as.vector(sapply(obs_bin_inds, function(bin_ind) {
      bin_ind + seq(from = -1, to = 1, by = 1)
    })))
    bins_to_sum <- bins_to_sum[
      bins_to_sum >= 1 & bins_to_sum <= length(expected_bin_names)]
  } else if(prediction_target %in% c("peak_inc", paste0("ph", 1:4, "_inc"))) {
    bins_to_sum <- obs_bin_inds + seq(from = -5, to = 5)
    bins_to_sum <- bins_to_sum[
      bins_to_sum >= 1 & bins_to_sum <= length(expected_bin_names)]
  }
  
  ## Do summation
  ## Futz around with bin names because order of expected bin names may not
  ## match order of bin_log_probs
  bin_names_to_sum <- expected_bin_names[bins_to_sum]
  log_prob <- logspace_sum(bin_log_probs[bin_names_to_sum])
  
  ## They truncate at -10
  log_prob <- max(-10, log_prob)
  
  ## return
  return(log_prob)
}



#' Get the onset baseline for a combination of region and season
#'
#' @param region a string, either "National", "Region k", or "Regionk" where
#'   k in {1, ..., 10}
#' @param season a string, in the format "2015/2016"
#'
#' @return baseline value for determining season onset
#'
#' @export
get_onset_baseline <- function(region, season = "2015/2016") {
  ## pick baseline
  ## assumes region is either "National" or "Region k" format
  reg_string <- ifelse(region=="National", "National", gsub(" ", "", region))
  idx <- which(cdcfluutils::flu_onset_baselines$region==reg_string &
                 cdcfluutils::flu_onset_baselines$season==season)
  reg_baseline <- cdcfluutils::flu_onset_baselines[idx, "baseline"]
  
  return(reg_baseline)
}



#' Given a set of sampled trajectories for future disease incidence, preprocess
#' them in preparation for using them to generate forecast samples for
#' short-term and seasonal incidence targets.
#' 
#' @param trajectory_samples an n by h matrix of sampled incidence trajectories;
#'     Each sampled trajectory is for the next h weeks. There are n trajectories
#' @param round_digits number of digits to which incidence will be rounded
#' @param obs_data a data frame containing observed incidence so far this
#'     season; at a minimum, must contain two columns: season_week, indicating
#'     the week of the current disease season, and the column given in
#'     prediction_target_var
#' @param prediction_target_var character specifying column name in obs_data for
#'     which we are generating predictions.
#' @param first_season_obs_ind index of row in obs_data for the first week of
#'     the season for which we are generating predictions.
#' @param analysis_time_ind index of row in obs_data that contains the most
#'     recent observation of incidence (default nrow(obs_data))
#' @param analysis_time_season in format 2018/2019
#' @param analysis_time_epiweek in format 201840
#' @param region
#' @param seasonal_target_week_limits vector of length 2 specifying season 
#'     weeks that will be used for identifying seasonal targets such as peak
#'    incidence; sampled or observed incidence outside of this range will be
#'    set to NA.
#' 
#' @return an n by W matrix of preprocessed trajectory samples, where
#'     W = (# of weeks observed so far this season) + h
#' 
#' @export
preprocess_and_augment_trajectory_samples <- function(
  trajectory_samples,
  round_digits,
  obs_data_matrix,
  obs_data,
  prediction_target_var,
  first_season_obs_ind,
  analysis_time_ind = nrow(obs_data),
  analysis_time_season,
  analysis_time_epiweek,
  region,
  seasonal_target_week_limits
) {
  ## Prepend observed data if supplied
  if(!missing(obs_data_matrix)) {
    trajectory_samples <- cbind(obs_data_matrix, trajectory_samples)
  }
  
  ## Round, if requested
  if(!missing(round_digits)) {
    trajectory_samples <- round(trajectory_samples, digits = round_digits)
  }
  
  ## If first observation for the season was not at season week 1,
  ## augment with leading NAs
  first_season_obs_week <- obs_data$season_week[first_season_obs_ind]
  if(first_season_obs_week != 1) {
    trajectory_samples <- cbind(
      matrix(NA, nrow = nrow(trajectory_samples), ncol = first_season_obs_week - 1),
      trajectory_samples
    )
  }
  
  ## values before the first analysis time week are NA so that
  ## onset and peak calculations only look at data within the CDC's definition
  ## of the flu season for purposes of the competition
  #trajectory_samples[, seq_len(seasonal_target_week_limits[1] - 1)] <- NA
  #trajectory_samples[,
  #  seasonal_target_week_limits[2] + seq_len(ncol(trajectory_samples) - seasonal_target_week_limits[2])] <- NA
  
  return(trajectory_samples)
}



#' Get samples in for each prediction target in a predx data frame
#'
#' for a given season using a predictive method that works by simulating
#' trajectories of incidence in each remaining week of the season.  Results are
#' stored in a data frame, saved in a .rds file with a name like
#' "model_name-region-season-loso-predictions.rds"
#' Results have columns indicating the analysis time season and season week,
#' model name, log scores for each prediction target, the "log score" used
#' in the competition (adding probabilities from adjacent bins) for each
#' prediction target, as well as the log of the probability assigned to each
#' bin.
#'
#' @param data data frame with observed disease data available to use for predictions
#' @param analysis_time_season character vector of length 1 specifying the
#'   season to obtain predictions for, in the format "2000/2001"
#' @param first_analysis_time_season_week integer specifying the first week of
#'   the season in which to make predictions, using all data up to and
#'   including that week to make predictions for each following week in the
#'   season
#' @param last_analysis_time_season_week integer specifying the last week of
#'   the season in which to make predictions, using all data up to and including
#'   that week to make predictions for each following week in the season
#' @param region string with name of region to make predictions for
#' @param prediction_target_var string specifying the name of the variable in
#'   data for which we want to make predictions
#' @param incidence_bins a data frame with variables lower and upper defining
#'   lower and upper endpoints to use in binning incidence
#' @param incidence_bin_names a character vector with a name for each incidence
#'   bin
#' @param n_trajectory_sims integer number of trajectories to simulate
#' @param simulate_trajectories_function a function to call to simulate
#'   incidence trajectories.  It will be called with the following arguments:
#'     * n_sims = number of trajectories to simulate
#'     * max_prediction_horizon = number of following weeks to simulate
#'     * data = all available data to use in doing simulation, up to and
#'         including the analysis time
#'     * region = region
#'     * analysis_time_season = analysis_time_season
#'     * analysis_time_season_week = week of the season at which we are making
#'         the predictions
#'     * params = simulate_trajectories_params; additional user-provided
#'         parameters
#' @param simulate_trajectories_params optional additional parameters to pass
#'   to simulate_trajectories_function
#' @param regional logical whether to make predictions for HHS regions (TRUE),
#'   or states (FALSE)
#'
#' @return predx data frame with samples for one region
#' @export
get_submission_one_region_via_trajectory_simulation <- function(
  data,
  analysis_time_season = "2019/2020",
  first_analysis_time_season_week = 10, # == week 40 of year
  last_analysis_time_season_week = 41, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
  region,
  prediction_target_var,
  incidence_bins,
  incidence_bin_names,
  n_trajectory_sims,
  simulate_trajectories_function,
  simulate_trajectories_params,
  backfill_adjust = c("none", "forecast_input", "post-hoc")) {
  
  backfill_adjust <- match.arg(backfill_adjust, choices = c("none", "forecast_input", "post-hoc"))
  
  weeks_in_first_season_year <-
    get_num_MMWR_weeks_in_first_season_year(analysis_time_season)
  
  ## subset data to be only the region-specific data
  data <- data[data$region == region,]
  
  analysis_time_season_week <- data$season_week[nrow(data)]
  analysis_time_ind <- nrow(data)
  max_prediction_horizon <- max(4L,
                                last_analysis_time_season_week + 1 - analysis_time_season_week)
  first_season_obs_ind <- min(which(data$season == analysis_time_season))
  
  if(backfill_adjust %in% c("post-hoc")) {
    epiweek <- cdcfluutils::season_week_to_year_week(
      analysis_time_season_week,
      weeks_in_first_season_year = weeks_in_first_season_year)
    std_region <- cdcfluutils::to_standard_location_code(region)
    obs_data_matrix <- rRevisedILI_cached(
      n = n_trajectory_sims,
      observed_inc = data[seq(from = first_season_obs_ind, to = analysis_time_ind), prediction_target_var, drop = TRUE],
      epiweek_idx = epiweek,
      region = std_region,
      season = analysis_time_season,
      season_start_epiweek = 31,
      add_nowcast = FALSE)
  } else if(backfill_adjust == "forecast_input") {
    epiweek <- cdcfluutils::season_week_to_year_week(
      analysis_time_season_week,
      weeks_in_first_season_year = weeks_in_first_season_year)
    std_region <- cdcfluutils::to_standard_location_code(region)
    temp <- rRevisedILI_cached(
      n = n_trajectory_sims,
      observed_inc = data[seq(from = first_season_obs_ind, to = analysis_time_ind), prediction_target_var, drop = TRUE],
      epiweek_idx = epiweek,
      region = std_region,
      season = analysis_time_season,
      season_start_epiweek = 31,
      add_nowcast = FALSE)
    obs_data_matrix <- temp$total_traj
    sampled_ids <- temp$sampled_id
  } else {
    obs_data_matrix <- matrix(
      rep(data[seq(from = first_season_obs_ind, to = analysis_time_ind), prediction_target_var, drop = TRUE],
          each = n_trajectory_sims),
      nrow = n_trajectory_sims
    )
  }
  
  if(identical(backfill_adjust, "forecast_input")) {
    trajectory_samples <- matrix(nrow = n_trajectory_sims, ncol = max_prediction_horizon)
    sampled_revised_data <- data
    unique_sampled_ids <- sampled_ids %>%
      count(region, season, epiweek)
    
    for(i in seq_len(nrow(unique_sampled_ids))) {
      inds <- which(
        sampled_ids$region == unique_sampled_ids$region[i] &
        sampled_ids$season == unique_sampled_ids$season[i] &
        sampled_ids$epiweek == unique_sampled_ids$epiweek[i]
      )
      sampled_revised_data[seq(from = first_season_obs_ind, to = analysis_time_ind), prediction_target_var] <- 
        obs_data_matrix[inds[1], ]
      trajectory_samples[inds, ] <- simulate_trajectories_function(
        n_sims = length(inds),
        max_prediction_horizon = max_prediction_horizon,
        data = sampled_revised_data[seq_len(analysis_time_ind), , drop = FALSE],
        region = region,
        analysis_time_season = analysis_time_season,
        analysis_time_season_week = analysis_time_season_week,
        params = simulate_trajectories_params
      )[seq_along(inds), ]
    }
  } else {
    trajectory_samples <- simulate_trajectories_function(
      n_sims = n_trajectory_sims,
      max_prediction_horizon = max_prediction_horizon,
      data = data[seq_len(analysis_time_ind), , drop = FALSE],
      region = region,
      analysis_time_season = analysis_time_season,
      analysis_time_season_week = analysis_time_season_week,
      params = simulate_trajectories_params
    )
  }
  
  # add observed incidence so far, round to nearest 0.1, and set incidence
  # outside CDC season bounds to NA
  trajectory_samples <- preprocess_and_augment_trajectory_samples(
    trajectory_samples = trajectory_samples,
    round_digits = 1,
    obs_data_matrix = obs_data_matrix,
    obs_data = data[seq_len(analysis_time_ind), , drop = FALSE],
    prediction_target_var = prediction_target_var,
    first_season_obs_ind = first_season_obs_ind,
    analysis_time_ind = analysis_time_ind,
    analysis_time_season = analysis_time_season,
    analysis_time_epiweek = cdcfluutils::season_week_to_year_week(analysis_time_season_week),
    region = region,
    seasonal_target_week_limits = c(first_analysis_time_season_week, last_analysis_time_season_week)
  )
  
  forecast_predx <- get_predx_forecasts_from_trajectory_samples(
    trajectory_samples = trajectory_samples,
    location = region,
    targets = c("Season onset", "Season peak week", "Season peak percentage",
      paste0(1:4, " wk ahead")),
    season = analysis_time_season,
    analysis_time_season_week = analysis_time_season_week,
    first_analysis_time_season_week = first_analysis_time_season_week,
    last_analysis_time_season_week = last_analysis_time_season_week,
    predx_types = c("Sample", "Bin", "Point")
  )
  
  return(forecast_predx)
}



#' @export
get_predx_forecasts_from_trajectory_samples <- function(
  trajectory_samples,
  location,
  targets = c("Season onset", "Season peak week", "Season peak percentage",
              paste0(1:4, " wk ahead")),
  season,
  analysis_time_season_week,
  first_analysis_time_season_week = 10, # == week 40 of year
  last_analysis_time_season_week = 41, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
  predx_types = c("Sample", "Bin", "Point")
) {
  result_targets <- NULL
  result_predx <- list()
  trajectory_samples <- as.matrix(trajectory_samples)
  
  if(any(c("Season onset", "Season peak week", "Season peak percentage") %in% targets)) {
    ## how many weeks in the first year of the season; used below
    weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(season)
    
    ## get indices in trajectory samples with NA values that affect
    ## estimation of seasonal quantities
    sample_inds_with_na <- apply(
      trajectory_samples[,
                         seq(from = first_analysis_time_season_week, to = last_analysis_time_season_week),
                         drop = FALSE],
      1,
      function(x) any(is.na(x)))
    
    ## Predictions for things about the whole season
    if(all(sample_inds_with_na)) {
      warning("NAs in all simulated trajectories, unable to predict seasonal quantities")
    } else {
      ## subset to sampled trajectories that are usable/do not have NAs
      subset_trajectory_samples <- trajectory_samples[!sample_inds_with_na, , drop = FALSE]
      
      ## values before the first analysis time week or after last are NA so
      ## that onset and peak calculations only look at data within the CDC's
      ## definition of the flu season for purposes of the competition
      subset_trajectory_samples[,
        seq(from = 1, to = first_analysis_time_season_week - 1)] <- NA_real_
      subset_trajectory_samples[,
        seq(from = last_analysis_time_season_week + 1,
          to = ncol(subset_trajectory_samples))] <- NA_real_
      
      ## Convert to binned values
      binned_subset_trajectory_samples <- get_inc_bin(
        subset_trajectory_samples, max = 13, return_character = FALSE)
      
      week_bins <- as.character(
        season_week_to_year_week(
          seq(from = 10, to = weeks_in_first_season_year - 10, by = 1),
          first_season_week = 31,
          weeks_in_first_season_year = weeks_in_first_season_year)
      )
      onset_week_bins <- c(week_bins, "none")
      
      ## Onset
      if("Season onset" %in% targets) {
        onset_week_by_sim_ind <-
          apply(binned_subset_trajectory_samples, 1, function(trajectory) {
            cdcfluutils::get_onset_week(
              incidence_trajectory = trajectory,
              baseline =
                cdcfluutils::get_onset_baseline(region = location, season = season),
              onset_length = 3L,
              first_season_week = 31,
              weeks_in_first_season_year = weeks_in_first_season_year
            )
          })
        onset_mmwr_week <- onset_week_by_sim_ind
        numeric_onset_inds <- which(onset_mmwr_week != "none")
        onset_mmwr_week[numeric_onset_inds] <- season_week_to_year_week(
              as.numeric(onset_week_by_sim_ind[numeric_onset_inds]),
              first_season_week = 31,
              weeks_in_first_season_year = weeks_in_first_season_year) %>%
            as.character()
        
        if("Sample" %in% predx_types) {
          result_targets <- c(result_targets, "Season onset")
          result_predx <- c(result_predx, predx::SampleCat(onset_mmwr_week))
        }
        if("Bin" %in% predx_types) {
          result_targets <- c(result_targets, "Season onset")
          result_predx <- c(result_predx,
            predx::transform_predx(predx::SampleCat(onset_mmwr_week),
              "BinCat",
              cat = onset_week_bins))
        }
        if("Point" %in% predx_types) {
          if(mean(onset_week_by_sim_ind == "none") < 1) {
            result_targets <- c(result_targets, "Season onset")
            pt_pred <- season_week_to_year_week(
                floor(median(suppressWarnings(as.numeric(onset_week_by_sim_ind)), na.rm = TRUE)),
                first_season_week = 31,
                weeks_in_first_season_year = weeks_in_first_season_year)
            result_predx <- c(result_predx, predx::Point(pt_pred))
          }
        }
      }
      
      ## Peak incidence
      if("Season peak percentage" %in% targets) {
        peak_inc_bin_by_sim_ind <-
          apply(binned_subset_trajectory_samples, 1, function(trajectory) {
            max(trajectory, na.rm = TRUE)
          })

        if("Sample" %in% predx_types) {
          result_targets <- c(result_targets, "Season peak percentage")
          result_predx <- c(result_predx, predx::Sample(peak_inc_bin_by_sim_ind))
        }
        if("Bin" %in% predx_types) {
          result_targets <- c(result_targets, "Season peak percentage")
          result_predx <- c(result_predx,
                            predx::transform_predx(predx::Sample(peak_inc_bin_by_sim_ind),
                                                   "BinLwr",
                                                   lwr = seq(from = 0.0, to = 13.0, by = 0.1)))
        }
        if("Point" %in% predx_types) {
          result_targets <- c(result_targets, "Season peak percentage")
          result_predx <- c(result_predx,
                            predx::Point(median(peak_inc_bin_by_sim_ind)))
        }
      }
      
      ## Peak week
      ## note that some sim inds may have more than 1 peak week...
      if("Season peak week" %in% targets) {
        # peak week is determined from rounded, but not binned, values
        peak_inc_bin_by_sim_ind <-
          apply(binned_subset_trajectory_samples, 1, function(trajectory) {
            max(trajectory, na.rm = TRUE)
          })
        
        peak_weeks_by_sim_ind <- unlist(lapply(
          seq_len(nrow(binned_subset_trajectory_samples)),
          function(sim_ind) {
            inc_val <- peak_inc_bin_by_sim_ind[sim_ind]
            peak_season_weeks <- which(
              binned_subset_trajectory_samples[sim_ind, ] == inc_val)
            return(peak_season_weeks)
          }
        ))
        
        peak_weeks_mmwr <- season_week_to_year_week(
          peak_weeks_by_sim_ind,
          first_season_week = 31,
          weeks_in_first_season_year = weeks_in_first_season_year)
        
        if("Sample" %in% predx_types) {
          result_targets <- c(result_targets, "Season peak week")
          result_predx <- c(result_predx, predx::SampleCat(as.character(peak_weeks_mmwr)))
        }
        if("Bin" %in% predx_types) {
          result_targets <- c(result_targets, "Season peak week")
          result_predx <- c(result_predx,
            predx::transform_predx(predx::SampleCat(as.character(peak_weeks_mmwr)),
              "BinCat",
              cat = week_bins))
        }
        if("Point" %in% predx_types) {
          result_targets <- c(result_targets, "Season peak week")
          result_predx <- c(result_predx,
            predx::Point(
              season_week_to_year_week(
                floor(median(peak_weeks_by_sim_ind)),
                first_season_week = 31,
                weeks_in_first_season_year = weeks_in_first_season_year)))
        }
      }
    }
  }
  
  ## Predictions for incidence in an individual week at prediction horizon ph = 1, ..., 4
  wk_targets <- targets[targets %in% paste0(1:4, " wk ahead")]
  phs <- as.numeric(substr(wk_targets, 1, 1))
  for(ph in phs) {
    if(analysis_time_season_week + ph <= ncol(trajectory_samples)) {
      sample_inds_with_na <- is.na(trajectory_samples[, analysis_time_season_week + ph])
      
      ## get sampled incidence values at prediction horizon that are usable/not NAs
      ph_inc_by_sim_ind <- trajectory_samples[!sample_inds_with_na, analysis_time_season_week + ph]
      ph_inc_bin_by_sim_ind <- get_inc_bin(ph_inc_by_sim_ind, return_character = FALSE)
      
      if("Sample" %in% predx_types) {
        result_targets <- c(result_targets, paste0(ph, " wk ahead"))
        result_predx <- c(result_predx, predx::Sample(ph_inc_bin_by_sim_ind))
      }
      if("Bin" %in% predx_types) {
        result_targets <- c(result_targets, paste0(ph, " wk ahead"))
        result_predx <- c(result_predx,
                          predx::transform_predx(predx::Sample(ph_inc_bin_by_sim_ind),
                                                 "BinLwr",
                                                 lwr = seq(from = 0.0, to = 13.0, by = 0.1)))
      }
      if("Point" %in% predx_types) {
        result_targets <- c(result_targets, paste0(ph, " wk ahead"))
        result_predx <- c(result_predx,
                          predx::Point(median(ph_inc_bin_by_sim_ind)))
      }
    }
  } # ph loop
  
  return(predx::as.predx_df(list(
    location = location,
    target = result_targets,
    predx = result_predx
  )))
}



#' return the bin name for a given incidence
#'
#' @param inc numeric incidence level
#' @param return_character logical: if true, return type is character (bin name)
#'   if false, return type is numeric representation of bin
#'
#' @return vector giving the bin name of the input incidence.
#'
#' @details assumes max inc bin is 13 and bins are 0.1 in size.
#'
#' @export
get_inc_bin <- function(inc,max=13,
                        return_character = TRUE) {
  inc <- round(inc, 1)
  bin_numeric <- ifelse(inc < max,
                        floor(inc*10)/10, ## floors to 1st decimal place
                        max)
  if(return_character) {
    return(as.character(bin_numeric))
  } else {
    return(bin_numeric)
  }
}



#' Calculation of median value from binned probability distribution
#'
#' @param probs vector of named probabilities
#'
#' @return a numeric value
#'
#' @export
calc_median_from_binned_probs <- function(probs) {
  ## could do something more intelligent for "none" bin in onset - currently assuming it is all ordered
  cumprob <- cumsum(probs)
  median_idx <- min(which(cumprob>=0.5))
  as.numeric(names(probs)[median_idx])
}

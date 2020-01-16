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
  reg_string <- cdcfluutils::to_standard_location_code(region)
  onset_baseline_reg_strings <- cdcfluutils::to_standard_location_code(
    cdcfluutils::flu_onset_baselines$region)
  
  idx <- which(
    onset_baseline_reg_strings == reg_string &
    cdcfluutils::flu_onset_baselines$season==season)
  
  reg_baseline <- cdcfluutils::flu_onset_baselines[idx, "baseline"]
  
  return(reg_baseline)
}



#' Find the first season week incidence goes below baseline and doesn't go back
#' above it at the end of the season
#' 
#' @parameter location character location code.  Will be converted to a
#'   standard location code via cdcfluutils::to_standard_location_code
#' @parameter season character season in format "2018/2019"
#' @parameter incidence_data data frame of incidence with at minimum columns
#'   `location`, `season`, `season_week`, and the variable specified in
#'   `target_variable``
#' @parameter target_variable character name of column with incidence
#'   (default "")
#' @parameter first_analysis_time_season_week the first season week for which we might score forecasts
#' @parameter last_analysis_time_season_week the last season week for which we might score forecasts
#' 
#' @return integer season week, or NA if incidence never went above baseline
#' 
#' @export
first_season_week_below_baseline <- function(
    location,
    season,
    incidence_data,
    target_variable,
    first_analysis_time_season_week = 10, # == epidemic week 40
    last_analysis_time_season_week =
      cdcfluutils::get_num_MMWR_weeks_in_first_season_year(season) - 11
) {
  onset_baseline <- cdcfluutils::get_onset_baseline(location, season)
  
  incidence_data_in_location_season <- incidence_data %>%
    filter(
      location == UQ(location),
      season == UQ(season),
      season_week >= first_analysis_time_season_week,
      season_week <= last_analysis_time_season_week
    )
  
  # after rounding to nearest 0.1, which indices are below baseline?
  below_baseline <-
    round(incidence_data_in_location_season[[target_variable]], 1) < onset_baseline
  
  if(!below_baseline[length(below_baseline)]) {
    # last observed incidence is at or above baseline
    return(last_analysis_time_season_week)
  } else if(all(below_baseline)) {
    # no observations above baseline
    return(NA_integer_)
  } else {
    # run length encoding of indices below/above baseline
    below_baseline_rle <- rle(
      incidence_data_in_location_season[[target_variable]] < onset_baseline
    )
    
    # keep all but last run (we know last run is below baseline)
    # the *following* week is the first below baseline
    inds_to_keep <- seq_len(length(below_baseline_rle$lengths) - 1)
    return(incidence_data_in_location_season$season_week[
      sum(below_baseline_rle$lengths[inds_to_keep])] + 1)
  }
}

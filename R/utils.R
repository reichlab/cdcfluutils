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
      ph_inc_bin_by_sim_ind <- get_inc_bin(ph_inc_by_sim_ind, return_character = FALSE) %>%
        as.numeric()
      
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

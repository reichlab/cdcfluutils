#' Preprocess state level flu data and set up a data object to pass to stan
#'
#' @param us_flu data frame of regional and national flu observations
#' @param seasons character vector of seasons to include in "2014/2015" format
#' @param predict_season which season are we generating predictions for
#' @param predict_last_obs_season_week integer season week from which to predict
#' @param return_extras boolean, default FALSE: return extra things not needed for
#'     stan, but that might be helpful for other things (like making plots)
#' @return list object with data for stan
#'
#' @export
preprocess_reg_data_for_GPthingy_stan <- function(
    us_flu,
    seasons,
    predict_season,
    predict_last_obs_season_week,
    return_extras = FALSE,
    knot_limits = c(1, 53L),
    DX = 10,
    get_all_region_seasons = FALSE,
    pred_all_t = FALSE) {
  # Subset to evaluation seasons, add variables with weeks
  us_flu <- us_flu %>%
    dplyr::filter(season %in% seasons) %>%
    dplyr::group_by(region) %>%
    dplyr::mutate(reg_week = dplyr::row_number())

  t0 <- us_flu$time[1]
  us_flu <- us_flu %>%
    mutate(weeks_since_t0 = as.numeric(time - t0) / 7 + 1)

  unique_regions <- unique(us_flu$region)
  S <- length(unique_regions)

  if(get_all_region_seasons) {
    unique_region_seasons <- expand.grid(region = unique_regions, season = seasons)
  } else {
    unique_region_seasons <- us_state_flu %>% dplyr::distinct(region, season)
  }

  SU <- nrow(unique_region_seasons)
  T_max <- 53L

  # Information about holidays
  holidays_considered <- c("thanksgiving", "christmas", "postchristmas")
  region_season_holiday_counts <- us_flu %>%
    dplyr::group_by(region, season) %>%
    dplyr::summarize_at(paste0(holidays_considered, "_week"), sum) %>%
    dplyr::ungroup()
  max_num_holidays <- max(region_season_holiday_counts %>% select(-region, -season))

  num_thanksgiving <-
    num_christmas <-
    num_postchristmas <-
      rep(max_num_holidays, SU)
  thanksgiving_week_inds <-
    christmas_week_inds <-
    postchristmas_week_inds <-
      rep(-1L, ncol = SU)

  y_t_combined <-
    n_t_combined <-
    ili_t_combined <-
    y_t_combined_full <-
    n_t_combined_full <-
    ili_t_combined_full <-
    a_t_combined <-
    p_t_combined <-
    prop_a_t_combined <-
    ili_obs_t_combined <-
    ili_obs_t_combined_full <-
    ili_pred_t_combined <-
    virology_obs_t_combined <-
      matrix(1.0,
        nrow = T_max,
        ncol = SU)

  num_obs_ili <- rep(0, SU)
  num_obs_ili_full <- rep(0, SU)
  num_pred_ili <- rep(0, SU)
  num_obs_virology <- rep(0, SU)

  for(su_ind in seq_len(SU)) {
    ## ili data for region season su
    s <- unique_region_seasons$region[su_ind]
    u <- unique_region_seasons$season[su_ind]
    region_season_data <- us_flu %>%
      dplyr::filter(region == s, season == u)

    if(u == predict_season && nrow(region_season_data) > 0) {
      t_present <- region_season_data$weeks_since_t0[
        region_season_data$season_week == predict_last_obs_season_week]
    } else {
      t_present <- Inf
    }

    y_t <- region_season_data$ilitotal
    n_t <- region_season_data$total_patients
    ili_t <- region_season_data$weighted_ili
    a_t <- region_season_data$virology_total_a
    p_t <- region_season_data$virology_total_a + region_season_data$virology_total_b
    prop_a_t <- a_t / p_t
    v_t <- region_season_data$virology_total_specimens
    obs_t <- as.numeric(region_season_data$season_week)
    weeks_since_t0 <- as.numeric(region_season_data$weeks_since_t0)

    ## drop missing values and observed 0's on holidays
    ili_non_na_inds <- (!is.na(y_t) & !is.na(n_t) & n_t > 0)
    #ili_non_na_inds <- ili_non_na_inds &
    #  !(y_t == 0 & region_season_data[[paste0(holidays_considered[1], "_week")]]) &
    #  !(y_t == 0 & region_season_data[[paste0(holidays_considered[2], "_week")]]) &
    #  !(y_t == 0 & region_season_data[[paste0(holidays_considered[3], "_week")]])

    y_t <- y_t[ili_non_na_inds]
    n_t <- n_t[ili_non_na_inds]
    ili_t <- ili_t[ili_non_na_inds]
    ili_obs_t <- obs_t[ili_non_na_inds]
    weeks_since_t0_ili <- weeks_since_t0[ili_non_na_inds]

    inds_to_t_present <- which(weeks_since_t0_ili <= t_present)
    y_t_to_t_present <- y_t[inds_to_t_present]
    n_t_to_t_present <- n_t[inds_to_t_present]
    ili_t_to_t_present <- ili_t[inds_to_t_present]
    ili_obs_t_to_t_present <- ili_obs_t[inds_to_t_present]

    y_t_combined[seq_along(y_t_to_t_present), su_ind] <- y_t_to_t_present
    n_t_combined[seq_along(n_t_to_t_present), su_ind] <- n_t_to_t_present
    ili_t_combined[seq_along(ili_t_to_t_present), su_ind] <- ili_t_to_t_present
    ili_obs_t_combined[seq_along(y_t_to_t_present), su_ind] <- ili_obs_t_to_t_present
    num_obs_ili[su_ind] <- length(y_t_to_t_present)

    y_t_combined_full[seq_along(y_t), su_ind] <- y_t
    n_t_combined_full[seq_along(n_t), su_ind] <- n_t
    ili_t_combined_full[seq_along(ili_t), su_ind] <- ili_t
    ili_obs_t_combined_full[seq_along(y_t), su_ind] <- ili_obs_t
    num_obs_ili_full[su_ind] <- length(y_t)

    if(u == predict_season) {
      ili_obs_t_to_predict <- ili_obs_t[-inds_to_t_present]
    } else {
      ili_obs_t_to_predict <- NULL
    }

    ili_pred_t_combined[seq_along(ili_obs_t_to_predict), su_ind] <- ili_obs_t_to_predict
    num_pred_ili[su_ind] <- length(ili_obs_t_to_predict)

    #thanksgiving_inds_used <- which(region_season_data$thanksgiving_week &
    #  ili_non_na_inds &
    #  region_season_data$weeks_since_t0 <= t_present)
    #thanksgiving_week_inds[seq_along(thanksgiving_inds_used), su_ind] <-
    #  thanksgiving_inds_used
    #num_thanksgiving[su_ind] <- length(thanksgiving_inds_used)

    christmas_week_inds[su_ind] <- as.numeric(region_season_data$season_week[region_season_data$christmas_week])
    #christmas_inds_used <- which(region_season_data$christmas_week &
    #  ili_non_na_inds &
    #  region_season_data$weeks_since_t0 <= t_present)
    #christmas_week_inds[seq_along(christmas_inds_used), su_ind] <-
    #  christmas_inds_used
    #num_christmas[su_ind] <- length(christmas_inds_used)

    #postchristmas_inds_used <- which(region_season_data$postchristmas_week &
    #  ili_non_na_inds &
    #  region_season_data$weeks_since_t0 <= t_present)
    #postchristmas_week_inds[seq_along(postchristmas_inds_used), su_ind] <-
    #  postchristmas_inds_used
    #num_postchristmas[su_ind] <- length(postchristmas_inds_used)

    virology_non_na_inds <- (!is.na(a_t) & !is.na(p_t) & p_t > 0)
    a_t <- a_t[virology_non_na_inds]
    p_t <- p_t[virology_non_na_inds]
    prop_a_t <- prop_a_t[virology_non_na_inds]
    virology_obs_t <- obs_t[virology_non_na_inds]
    weeks_since_t0_virology <- weeks_since_t0[virology_non_na_inds]

    inds_to_t_present <- which(weeks_since_t0_virology <= t_present)
    a_t <- a_t[inds_to_t_present]
    p_t <- p_t[inds_to_t_present]
    prop_a_t <- prop_a_t[inds_to_t_present]
    virology_obs_t <- virology_obs_t[inds_to_t_present]

    a_t_combined[seq_along(a_t), su_ind] <- a_t
    p_t_combined[seq_along(a_t), su_ind] <- p_t
    prop_a_t_combined[seq_along(a_t), su_ind] <- prop_a_t
    virology_obs_t_combined[seq_along(a_t), su_ind] <- virology_obs_t
    num_obs_virology[su_ind] <- length(a_t)
  }

  all_knots <- seq(from = knot_limits[1], to = knot_limits[2], length = DX)
  center_knots <- all_knots[seq(from = 2, to = DX - 1)]
  X <- splines::bs(
    x = seq(from = 1, to = T_max, by = 1),
    knots = center_knots,
    intercept = TRUE
  )
  X_for_basis_centers <- splines::bs(
    x = seq(from = knot_limits[1], to = knot_limits[2], by = 1),
    knots = center_knots,
    intercept = TRUE
  )

  ili_t_concat <- unlist(lapply(seq_len(SU),
  function(su_ind) {
    y_t_combined[seq_len(num_obs_ili[su_ind]), su_ind] / n_t_combined[seq_len(num_obs_ili[su_ind]), su_ind]
  }))

  ili_t_full_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      y_t_combined_full[seq_len(num_obs_ili_full[su_ind]), su_ind] / n_t_combined_full[seq_len(num_obs_ili_full[su_ind]), su_ind]
    }))

  ili_obs_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      ili_obs_t_combined[seq_len(num_obs_ili[su_ind]), su_ind]
    }))

  ili_obs_t_full_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      ili_obs_t_combined_full[seq_len(num_obs_ili_full[su_ind]), su_ind]
    }))

  ili_pred_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      ili_pred_t_combined[seq_len(num_pred_ili[su_ind]), su_ind]
    }))

  y_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      y_t_combined[seq_len(num_obs_ili[su_ind]), su_ind]
    }))

  n_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      n_t_combined[seq_len(num_obs_ili[su_ind]), su_ind]
    }))

  ili_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      ili_t_combined[seq_len(num_obs_ili[su_ind]), su_ind]
    }))

  virology_obs_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      virology_obs_t_combined[seq_len(num_obs_virology[su_ind]), su_ind]
    }))

  a_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      a_t_combined[seq_len(num_obs_virology[su_ind]), su_ind]
    }))

  p_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      p_t_combined[seq_len(num_obs_virology[su_ind]), su_ind]
    }))

  prop_a_t_concat <- unlist(lapply(seq_len(SU),
    function(su_ind) {
      prop_a_t_combined[seq_len(num_obs_virology[su_ind]), su_ind]
    }))

  obs_inds <- c()
  pred_inds <- c()
  for(su in seq_len(SU)) {
    base_ind <- sum(num_obs_ili[seq_len(su - 1)])
    base_ind_full <- sum(num_obs_ili_full[seq_len(su - 1)])
    obs_inds <- c(obs_inds,
      base_ind_full + which(
        ili_obs_t_full_concat[seq(from = base_ind_full + 1, to = base_ind_full + num_obs_ili_full[su])] %in%
          ili_obs_t_concat[seq(from = base_ind + 1, to = base_ind + num_obs_ili[su])]
      )
    )
    pred_inds <- c(pred_inds,
      base_ind_full + which(!(
        ili_obs_t_full_concat[seq(from = base_ind_full + 1, to = base_ind_full + num_obs_ili_full[su])] %in%
          ili_obs_t_concat[seq(from = base_ind + 1, to = base_ind + num_obs_ili[su])]
      ))
    )
  }

  unique_regions <- unique(unique_region_seasons$region)
  unique_seasons <- unique(unique_region_seasons$season)
  S <- length(unique_regions)
  U_max <- length(unique_seasons)

  SU_to_S_map <- sapply(
    seq_len(SU),
    function(su) {which(unique_regions == unique_region_seasons$region[su])}
  )

  SU_to_U_map <- sapply(
    seq_len(SU),
    function(su) {which(unique_seasons == unique_region_seasons$season[su])}
  )

  stan_data <- list(
    SU = SU,
    S = S,
    U_max = U_max,
    SU_to_S_map = SU_to_S_map,
    SU_to_U_map = SU_to_U_map,
    T_max = T_max,
    #max_num_holidays = max_num_holidays,
    #num_thanksgiving = sum(num_thanksgiving),
    #num_christmas = sum(num_christmas),
    #num_postchristmas = sum(num_postchristmas),
    #num_thanksgiving_by_state = num_thanksgiving,
    #num_christmas_by_state = num_christmas,
    #num_postchristmas_by_state = num_postchristmas,
    #thanksgiving_week_inds = thanksgiving_week_inds,
    christmas_week_ind = christmas_week_inds,
    #postchristmas_week_inds = postchristmas_week_inds,
    num_obs_ili = num_obs_ili,
    total_num_obs_ili = sum(num_obs_ili),
    #y_t = y_t_combined,
    #n_t = n_t_combined,
    y_t_concat = y_t_concat[seq_len(sum(num_obs_ili))],
    n_t_concat = n_t_concat[seq_len(sum(num_obs_ili))],
    ili_t_concat = ili_t_concat[seq_len(sum(num_obs_ili))] / 100.0,
    #ili_t_concat = ili_t_concat[seq_len(sum(num_obs_ili[SU_inds_to_use]))] * 100 + 0.05, #y_t_combined / n_t_combined,
    ili_obs_t = ili_obs_t_concat[seq_len(sum(num_obs_ili))],
    num_pred_ili = num_pred_ili,
    total_num_pred_ili = sum(num_pred_ili),
    ili_pred_t = ili_pred_t_concat[seq_len(sum(num_pred_ili))],
    num_obs_virology = num_obs_virology,
    total_num_obs_virology = sum(num_obs_virology),
    a_t_concat = a_t_concat[seq_len(sum(num_obs_virology))],
    p_t_concat = p_t_concat[seq_len(sum(num_obs_virology))],
    prop_a_t_concat = prop_a_t_concat[seq_len(sum(num_obs_virology))],
    virology_obs_t = virology_obs_t_concat[seq_len(sum(num_obs_virology))],
    #num_obs_virology = num_obs_virology,
    #a_t = a_t_combined,
    #p_t = p_t_combined,
    #virology_obs_t = virology_obs_t_combined,
    DX = ncol(X),
    X = X,
    knots = apply(X_for_basis_centers, 2, function(x) {weighted.mean(seq(from = knot_limits[1], to = knot_limits[2]), w = x)}),
    #knots = all_knots,
    max_lag = 1,
    gamma_trans = 0.0,
    lambda_trans = 0)#car::powerTransform(matrix(ili_t_concat[seq_len(sum(num_obs_ili[SU_inds_to_use]))] * 100 + 0.5), family = "bcPower")$lambda)

  if(return_extras) {
    stan_data$num_obs_ili_full = num_obs_ili_full
  }

  return(stan_data)
}

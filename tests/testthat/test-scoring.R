context("scoring")
library(dplyr)
library(cdcfluutils)


test_that("first_season_week_below_baseline correct", {
  # all below baseline returns NA_integer_
  location <- "hhs region 1"
  season <- "2017/2018"
  test_input <- data.frame(
    location = location,
    season = season,
    season_week = 1:52,
    weighted_ili = rep(cdcfluutils::get_onset_baseline(location, season) - 0.1, 52)
  )
  
  expected <- NA_integer_
  actual <- first_season_week_below_baseline(
    location = location,
    season = season,
    incidence_data = test_input,
    target_variable = "weighted_ili"
  )
  expect_identical(expected, actual)
  
  # all at or above baseline returns last season week scored
  location <- "hhs region 1"
  season <- "2017/2018"
  test_input <- data.frame(
    location = location,
    season = season,
    season_week = 1:52,
    weighted_ili = rep(cdcfluutils::get_onset_baseline(location, season), 52)
  )
  
  expected <- 41.0
  actual <- first_season_week_below_baseline(
    location = location,
    season = season,
    incidence_data = test_input,
    target_variable = "weighted_ili"
  )
  expect_identical(expected, actual)
  
  # some at or above baseline and some below returns first season week below
  # baseline in last run below baseline
  location <- "hhs region 1"
  season <- "2017/2018"
  obs_weighted_ili <- rep(cdcfluutils::get_onset_baseline(location, season), 52)
  obs_weighted_ili[36:52] <- obs_weighted_ili[36:52] - 0.1
  obs_weighted_ili[25:30] <- obs_weighted_ili[25:30] - 0.1
  
  test_input <- data.frame(
    location = location,
    season = season,
    season_week = 1:52,
    weighted_ili = obs_weighted_ili
  )
  
  expected <- 36.0
  actual <- first_season_week_below_baseline(
    location = location,
    season = season,
    incidence_data = test_input,
    target_variable = "weighted_ili"
  )
  expect_identical(expected, actual)
})

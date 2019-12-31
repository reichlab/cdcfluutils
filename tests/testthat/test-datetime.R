context("datetime")
library(dplyr)
library(cdcfluutils)


test_that("mmwr_week_to_season correct", {
  # separate years provided
  mmwr_weeks <- c(40, 30, 31, 20)
  mmwr_years <- c(2017, 2018, 2019, 2020)
  expected_seasons <- c("2017/2018", "2017/2018", "2019/2020", "2019/2020")
  
  expect_identical(
    cdcfluutils::mmwr_week_to_season(mmwr_weeks, mmwr_years, first_season_week = 31),
    expected_seasons
  )
  
  # one year provided
  mmwr_weeks <- c(40, 30, 31, 20)
  mmwr_years <- c(2017)
  expected_seasons <- c("2017/2018", "2016/2017", "2017/2018", "2016/2017")
  
  expect_identical(
    cdcfluutils::mmwr_week_to_season(mmwr_weeks, mmwr_years, first_season_week = 31),
    expected_seasons
  )
})

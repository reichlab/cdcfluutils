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



test_that("season_week_to_mmwr_year correct", {
  # first year, 52 week season
  actual <- cdcfluutils::season_week_to_mmwr_year(
    season_week = 21,
    season = "2018/2019",
    first_season_week = 31)
  expected <- 2018L
  expect_identical(expected, actual)

  # second year, 52 week season
  actual <- cdcfluutils::season_week_to_mmwr_year(
    season_week = 22,
    season = "2018/2019",
    first_season_week = 31)
  expected <- 2019L
  expect_identical(expected, actual)

  # first year, 53 week season
  actual <- cdcfluutils::season_week_to_mmwr_year(
    season_week = 21,
    season = "2014/2015",
    first_season_week = 31)
  expected <- 2014L
  expect_identical(expected, actual)

  # first year, 53 week season
  actual <- cdcfluutils::season_week_to_mmwr_year(
    season_week = 22,
    season = "2014/2015",
    first_season_week = 31)
  expected <- 2014L
  expect_identical(expected, actual)

  # first year, 53 week season
  actual <- cdcfluutils::season_week_to_mmwr_year(
    season_week = 23,
    season = "2014/2015",
    first_season_week = 31)
  expected <- 2015L
  expect_identical(expected, actual)
})


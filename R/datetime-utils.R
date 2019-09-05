#' Convert season week to year week
#'
#' @param season_week vector of indices of weeks in season
#' @param first_season_week number of week in year corresponding to the first
#'   week in the season.  For example, our code takes this value to be 31:
#'   a new influenza season starts on the 31st week of each year.
#' @param weeks_in_first_season_year How many MMWR weeks are in the first year
#'   of the season?  For example, in the 2000/2001 season, the first year is
#'   2000.  There were 52 MMWR weeks in 2000.
#'
#' @return vector of the same length as season_week with the week of the year
#'   that each observation falls on
#'
#' @export
season_week_to_year_week <- function(
  season_week,
  first_season_week = 31,
  weeks_in_first_season_year) {
  
  year_week <- season_week
  
  ## For competition first bin is week 40
  year_week[season_week < 10] <- 40
  year_week[season_week >= 10] <- season_week + first_season_week - 1
  year_week[year_week > weeks_in_first_season_year] <-
    year_week[year_week > weeks_in_first_season_year] -
    weeks_in_first_season_year
  
  return(year_week)
}



#' Convert year week to season week
#'
#' @param year_week vector of indices of weeks in year
#' @param year either a single (four digit) year or a vector of years with the
#'   same length as year_week
#'
#' @return vector of the same length as year_week with the week of the season
#'   that each observation falls on
#'
#' @export
year_week_to_season_week <- function(
  year_week,
  year) {
  season_week <- ifelse(
    year_week <= 30,
    year_week + MMWRweek::MMWRweek(MMWRweek:::start_date(year) - 1)$MMWRweek - 30,
    year_week - 30
  )
  
  return(season_week)
}



#' return integer that's either 52 or 53: number of MMWR weeks in the given year
#'
#' @param year year in the format "2014" -- can be character or numeric
#'
#' @details requires non-exported function start_date from MMWRweek package
#'
#' @export
get_num_MMWR_weeks_in_year <- function(year) {
  require(MMWRweek)
  year <- as.numeric(year)
  return(MMWRweek::MMWRweek(MMWRweek:::start_date(year + 1) - 1)$MMWRweek)
}



#' return integer that's either 52 or 53: number of weeks in the first year of
#' a given season.
#'
#' @param season season in the format "2014/2015"
#'
#' @details requires MMWRweek package
#'
#' @export
get_num_MMWR_weeks_in_first_season_year <- function(season) {
  return(get_num_MMWR_weeks_in_year(substr(season, 1, 4)))
}

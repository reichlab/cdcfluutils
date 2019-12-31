
#' For each week specified by a vector of dates giving the week start dates,
#' determine whether the week contains a specified date.
#'
#' @param week_start_date A vector of Date objects specifying the date of the first day in
#' the weeks of interest
#' @param year_to_pick integer or character giving the year to pick, e.g. "2010"
#' @param month_to_pick integer or character giving the month to pick, e.g. "12"
#' @param day_to_pick integer or character giving the day to pick, e.g. "22"
#'
#' @return a logical vector of the same length as time.  Entry i is TRUE if the
#' week beginning on week_start_date[i] contains the date specified by
#' year_to_pick, month_to_pick, and day_to_pick; FALSE otherwise.
#'
#' @export
pick_week <- function(
  week_start_date,
  year_to_pick,
  month_to_pick,
  day_to_pick) {
  date_to_pick <- lubridate::ymd(paste(year_to_pick, month_to_pick, day_to_pick, sep = "-"))
  selTF <- (week_start_date >= date_to_pick - 6) & (week_start_date <= date_to_pick)
  return(selTF)
}


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



#' get season in which an mmwr week falls
#' 
#' @param mmwr_week mmwr week as an integer (only the week)
#' @param mmwr_year year for mmwr week as an integer
#' first_season_week = 31,
#' @param first_season_week number of week in year corresponding to the first
#'   week in the season.  For example, our code takes this value to be 31:
#'   a new influenza season starts on the 31st week of each year.
#' 
#' @return vector of the same length as mmwr_week with season in which the week
#' falls, in the format "2018/2019"
#' 
#' @export
mmwr_week_to_season <- function(
  mmwr_week,
  mmwr_year,
  first_season_week = 31) {
  season <- ifelse(
    mmwr_week < first_season_week,
    paste0(mmwr_year - 1, "/", mmwr_year),
    paste0(mmwr_year, "/", mmwr_year + 1)
  )
  
  return(season)
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

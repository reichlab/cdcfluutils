#' Wrapper function around epiforecast::mimicPastEpidataDF to return partially
#' revised CDC flu data available as of the specified epiweek (but with missing
#' values filled in with more recent data from ilinet).  On top of
#' epiforecast::mimicPastEpidataDF, this function: 1) uses a static pull of
#' data from the delphi API, stored in the data directory of this package;
#' 2) subsets data to a specified region; and 3) drops some occasional extra
#' rows of NA's at the end of the output from epiforecast::mimicPastEpidataDF.
#'
#' @param region_str character defining region of interest, must be in c("nat", paste0("hhs", 1:10))
#' @param epiweek_str character string defining an epiweek in YYYYWW format
#'
#' @return a dataset in similar format to that returned by the Delphi epidata API
#'
#'
#' @export
get_partially_revised_ilinet <- function(region_str, epiweek_str) {
  require(epiforecast)
  require(dplyr)

  ## from ilinet, via DELPHI API
  partially_revised_data <- readRDS(file.path(
      find.package("cdcfluutils"),
      "data",
      "flu_data_with_backfill.rds")) %>%
    dplyr::filter(region == region_str)

  temp <- epiforecast::mimicPastEpidataDF(
      partially_revised_data,
      epiweek_str) %>%
    dplyr::filter(epiweek <= as.integer(epiweek_str)) %>%
    as.data.frame()

  return(temp)
}


#' Download backfill data from the delphi api
#' and save to data folder.
#' Note: We may want to do this periodically to make sure we are up to date
download_backfill_data <- function(){
  library(plyr) # for rbind.fill
  library(dplyr)
  source("https://raw.githubusercontent.com/cmu-delphi/delphi-epidata/master/src/client/delphi_epidata.R")
  
  # Fetch data
  all_obs <- lapply(c("nat", paste0("hhs", 1:10)),
                    function(region_val) {
                      lapply(1:51,
                             function(lag_val) {
                               obs_one_lag <- Epidata$fluview(
                                 regions = list(region_val),
                                 epiweeks = list(Epidata$range(199740, 201815)),
                                 lag = list(lag_val))
                               
                               lapply(obs_one_lag$epidata,
                                      function(x) {
                                        x[sapply(x, function(comp) is.null(comp))] <- NA
                                        return(as.data.frame(x))
                                      }) %>%
                                 rbind.fill()
                             }) %>%
                        rbind.fill()
                    }) %>%
    rbind.fill()
  
  saveRDS(all_obs,
          file = "data/flu_data_with_backfill.rds")
  
}



#' Download and preprocess the latest CDC flu data, both national and regional
#'
#' @param latest_year year through which data should be downloaded, defaults to current year
#'
#' @return data frame with latest flu data, preprocessed
#' @export
download_and_preprocess_flu_data <- function(latest_year = as.numeric(format(Sys.Date(), "%Y"))) {
  require(cdcfluview)
  require(lubridate)
  require(dplyr)
  require(MMWRweek)
  
  regionflu <- ilinet(region="hhs", years= 1997:latest_year)
  regionflu$region <- as.character(regionflu$region)
  
  usflu <- ilinet(region="national", years= 1997:latest_year)
  
  flu_data <- bind_rows(regionflu, usflu)
  
  flu_data <- transmute(flu_data,
                        region_type = region_type,
                        region = as.factor(region),
                        year = year,
                        week = week,
                        time = as.POSIXct(MMWRweek2Date(year, week)),
                        weighted_ili = weighted_ili)
  
  ## set zeroes to NAs
  flu_data[which(flu_data$weighted_ili==0),"weighted_ili"] <- NA
  
  ## Add time_index column: the number of days since some origin date
  ## (1970-1-1 in this case).  The origin is arbitrary.
  flu_data$time_index <- as.integer(lubridate::date(flu_data$time) -  ymd("1970-01-01"))
  
  ## Season column: for example, weeks of 2010 up through and including week 30
  ## get season 2009/2010; weeks after week 30 get season 2010/2011
  ## Official CDC flu season for the purposes of prediction runs from week 40 of
  ## one year to week 20 of the next; the season start week we define here is the
  ## mid-point of the "off-season"
  flu_data$season <- ifelse(
    flu_data$week <= 30,
    paste0(flu_data$year - 1, "/", flu_data$year),
    paste0(flu_data$year, "/", flu_data$year + 1)
  )
  
  ## Season week column: week number within season
  ## weeks after week 30 get season_week = week - 30
  ## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
  ## This computation relies on the start_date function in package MMWRweek,
  ## which is not exported from that package's namespace!!!
  flu_data$season_week <- ifelse(
    flu_data$week <= 30,
    flu_data$week + MMWRweek(MMWRweek:::start_date(flu_data$year) - 1)$MMWRweek - 30,
    flu_data$week - 30
  )
  
  flu_data <- as.data.frame(flu_data)
  
  return(flu_data)
}



#' Download and preprocess the latest CDC flu data, state-level
#'
#' @param latest_year year through which data should be downloaded, defaults to current year
#'
#' @return data frame with latest state-level flu data, preprocessed
#' @export
download_and_preprocess_state_flu_data <- function(latest_year = as.numeric(format(Sys.Date(), "%Y"))) {
  
  require(cdcfluview)
  require(MMWRweek)
  require(dplyr)
  require(lubridate)
  
  flu_data_raw <- ilinet(region="state", years=1997:latest_year)
  
  flu_data <- mutate(flu_data_raw, time = as.POSIXct(MMWRweek2Date(year, week)))
  
  ## set rows with denominator zeroes to NAs
  flu_data[which(flu_data$total_patients==0),"weighted_ili"] <- NA
  
  ## Add time_index column: the number of days since some origin date
  ## (1970-1-1 in this case).  The origin is arbitrary.
  flu_data$time_index <- as.integer(lubridate::date(flu_data$time) -  ymd("1970-01-01"))
  
  ## Season column: for example, weeks of 2010 up through and including week 30
  ## get season 2009/2010; weeks after week 30 get season 2010/2011
  ## Official CDC flu season for the purposes of prediction runs from week 40 of
  ## one year to week 20 of the next; the season start week we define here is the
  ## mid-point of the "off-season"
  flu_data$season <- ifelse(
    flu_data$week <= 30,
    paste0(flu_data$year - 1, "/", flu_data$year),
    paste0(flu_data$year, "/", flu_data$year + 1)
  )
  
  ## Season week column: week number within season
  ## weeks after week 30 get season_week = week - 30
  ## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
  ## This computation relies on the start_date function in package MMWRweek,
  ## which is not exported from that package's namespace!!!
  flu_data$season_week <- ifelse(
    flu_data$week <= 30,
    flu_data$week + MMWRweek(MMWRweek:::start_date(flu_data$year) - 1)$MMWRweek - 30,
    flu_data$week - 30
  )
  
  state_flu <- as.data.frame(flu_data)
  
  return(state_flu)
}



#' Download and preprocess the latest CDC hospitalization data
#'
#' @param latest_year year through which data should be downloaded, defaults to current year
#'
#' @return data frame with latest hospitalization data, preprocessed
#' @export
download_and_preprocess_hosp_data <- function(latest_year = as.numeric(format(Sys.Date(), "%Y"))) {
  
  require(cdcfluview)
  require(MMWRweek)
  require(dplyr)
  require(lubridate)
  
  flu_data_raw_hosp <-cdcfluview::hospitalizations(years=2009:latest_year)
  flu_data_raw_hosp$week <- flu_data_raw_hosp$year_wk_num
  flu_data <- mutate(flu_data_raw_hosp, time = as.POSIXct(MMWRweek2Date(year, week)))
  
  ## set rows with denominator zeroes to NAs
  #flu_data[which(flu_data$total_patients==0),"weighted_ili"] <- NA
  
  ## Add time_index column: the number of days since some origin date
  ## (1970-1-1 in this case).  The origin is arbitrary.
  flu_data$time_index <- as.integer(lubridate::date(flu_data$time) -  ymd("1970-01-01"))
  
  ## Season column: for example, weeks of 2010 up through and including week 30
  ## get season 2009/2010; weeks after week 30 get season 2010/2011
  ## Official CDC flu season for the purposes of prediction runs from week 40 of
  ## one year to week 20 of the next; the season start week we define here is the
  ## mid-point of the "off-season"
  flu_data$season <- ifelse(
    flu_data$week <= 30,
    paste0(flu_data$year - 1, "/", flu_data$year),
    paste0(flu_data$year, "/", flu_data$year + 1)
  )
  
  ## Season week column: week number within season
  ## weeks after week 30 get season_week = week - 30
  ## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
  ## This computation relies on the start_date function in package MMWRweek,
  ## which is not exported from that package's namespace!!!
  flu_data$season_week <- ifelse(
    flu_data$week <= 30,
    flu_data$week + MMWRweek(MMWRweek:::start_date(flu_data$year) - 1)$MMWRweek - 30,
    flu_data$week - 30
  )
  
  hosp_flu <- as.data.frame(flu_data)
  
  return(hosp_flu)
}

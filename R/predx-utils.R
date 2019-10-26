#' merge predx samples by location and target
#' 
#' @param predx_df a predx data frame with multiple rows that have samples at
#'   the same location and target
#' @return a single predx data frame with merged samples.  Only rows with
#'   predx_class Sample or SampleCat are retained.
#'   
#' @export
merge_predx_samples <- function(predx_df) {
  unique_location_targets <- predx_df %>%
    distinct(location, target)
  
  merged_predx <- purrr::pmap_dfr(
    unique_location_targets,
    function(location, target) {
      this_location <- location
      this_target <- target
      predx_samples <- predx_df %>%
        filter(predx_class %in% c("Sample", "SampleCat"),
               location == this_location, target == this_target)
      
      if(predx_samples$predx_class[1] == "Sample") {
        all_sample_vals <- unlist(lapply(predx_samples$predx, function(x) {x@predx}))
        predx_sample <- new("Sample", predx = all_sample_vals)
      } else {
        all_sample_vals <- unlist(lapply(predx_samples$predx, function(x) {x@predx}))
        predx_sample <- new("SampleCat", predx = all_sample_vals)
      }
      
      result <- predx::as.predx_df(list(
        location = this_location,
        target = this_target,
        predx = list(predx_sample)
      ))
      
      return(result)
    }
  )
  
  return(merged_predx)
}


#' Converts predx format forecasts into a data frame suitable for writing out
#' a csv file for submission
#' 
#' @param predx_df predx data frame object with samples for submission
#' @param year year, numeric or character
#' @param ew epidemic week, character in format "01"
#' @param team team name, character
#' 
#' @return path to csv file
#' 
#' @export
predx_to_submission_df <- function(predx_df, ew, year, team = "Kernel of Truth") {
  year <- as.numeric(year)
  if(as.numeric(ew) <= 30) {
    season <- paste0(year - 1, "/", year)
  } else {
    season <- paste0(year, "/", year + 1)
  }
  weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(season)
  
  unique_location_targets <- predx_df %>%
    distinct(location, target)
  
  augmented_predx <- purrr::pmap_dfr(
    unique_location_targets,
    function(location, target) {
      temp <- predx_df %>%
        filter(predx_class %in% c("Sample", "SampleCat"),
          location == UQ(location), target == UQ(target))
      
      if(nrow(temp) != 1) {
        stop("predx_df should have exactly 1 Sample or SampleCat row for each combination of location and target")
      }
      
      if(temp$predx_class == "Sample") {
        predx_sample <- temp$predx[[1]]
        
        predx_bin <- predx::transform_predx(predx_sample,
          "BinLwr",
          lwr = seq(from = 0.0, to = 13.0, by = 0.1))
      } else {
        predx_sample <- temp$predx[[1]]
        
        week_bins <- as.character(
          season_week_to_year_week(
            seq(from = 10, to = weeks_in_first_season_year - 10, by = 1),
            first_season_week = 31,
            weeks_in_first_season_year = weeks_in_first_season_year)
        )
        onset_week_bins <- c(week_bins, "none")
        
        if(target == "Season onset") {
          predx_bin <- predx::transform_predx(
            predx_sample,
            "BinCat",
            cat = onset_week_bins)
        } else if(target == "Season peak week") {
          predx_bin <- predx::transform_predx(
            predx_sample,
            "BinCat",
            cat = week_bins
          )
        }
      }
      
      result <- predx::as.predx_df(list(
        location = rep(location, 2),
        target = rep(target, 2),
        predx = list(predx_sample, predx_bin)
      ))

      return(result)
    }
  )

  submission_df <- augmented_predx %>%
   dplyr::filter(predx_class %in% c("BinCat", "BinLwr", "Point", "PointCat")) %>%
   dplyr::mutate(
     team = team,
     mmwr_week = paste0(year, ew),
     submission_date = Sys.Date(),
     unit = ifelse(target %in% c("Season onset", "Season peak week"), "week", "percent")
   ) %>%
   predx::export_flusight_csv()
  
  point_forecasts <- submission_df %>%
   FluSight::generate_point_forecasts()
  
  submission_df <- bind_rows(
    submission_df,
    point_forecasts) %>%
    arrange(location, target, desc(type))
  
  return(submission_df)
}
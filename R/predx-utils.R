#' Combines predx format forecast files from multiple different runs and writes
#' the results to csv file.
#' 
#' @param file_path file path to base directory for method.  should contain
#'  subfolder predx with .rds files with predx format files and subfolder csv
#'  into which results will be saved
#' @param ew epidemic week, character in format "01"
#' @param year year, numeric or character
#' 
#' @return path to csv file
#' 
#' @export
predx_to_csv <- function(file_path, ew, year) {
  files <- Sys.glob(paste0(file_path, "predx/EW", ew, "-", year, "*.rds"))
  year <- as.numeric(year)
  
  if(length(files) == 1) {
    merged_predx <- readRDS(files)
    files <- gsub("post-hoc", "post_hoc", files)
    method <- tail(strsplit(substr(files, 1, nchar(files) - 4), "-")[[1]], 1)
  } else {
    filename_for_method <- gsub("post-hoc", "post_hoc", files[1])
    method <- tail(strsplit(substr(filename_for_method, 1, nchar(filename_for_method) - 4), "-")[[1]], 2)[[1]]
    if(as.numeric(ew) <= 30) {
      season <- paste0(year - 1, "/", year)
    } else {
      season <- paste0(year, "/", year + 1)
    }
    weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(season)
    
    all_predx <- purrr::map_dfr(files, readRDS)
    
    unique_location_targets <- all_predx %>%
      distinct(location, target)
    
    merged_predx <- purrr::pmap_dfr(
      unique_location_targets,
      function(location, target) {
        this_location <- location
        this_target <- target
        temp <- all_predx %>%
          filter(predx_class %in% c("Sample", "SampleCat"),
            location == this_location, target == this_target)
        
        if(temp$predx_class[1] == "Sample") {
          all_sample_vals <- unlist(lapply(temp$predx, function(x) {x@predx}))
          predx_sample <- new("Sample", predx = all_sample_vals)
          
          predx_bin <- predx::transform_predx(predx::Sample(as.numeric(all_sample_vals)),
            "BinLwr",
            lwr = seq(from = 0.0, to = 13.0, by = 0.1))
          
          predx_point <- predx::Point(median(as.numeric(all_sample_vals)))
        } else {
          all_sample_vals <- unlist(lapply(temp$predx, function(x) {x@predx}))
          predx_sample <- new("SampleCat", predx = all_sample_vals)
          
          week_bins <- as.character(
            season_week_to_year_week(
              seq(from = 10, to = weeks_in_first_season_year - 10, by = 1),
              first_season_week = 31,
              weeks_in_first_season_year = weeks_in_first_season_year)
          )
          onset_week_bins <- c(week_bins, "none")
          
          if(this_target == "Season onset") {
            predx_bin <- predx::transform_predx(
              predx::SampleCat(all_sample_vals),
              "BinCat",
              cat = onset_week_bins)
            predx_point <- predx::PointCat("NA")
          } else if(this_target == "Season peak week") {
            predx_bin <- predx::transform_predx(
              predx::SampleCat(all_sample_vals),
              "BinCat",
              cat = week_bins
            )
            predx_point <- predx::PointCat("NA")
          }
        }
        
        if(is.null(predx_point)) {
          result <- predx::as.predx_df(list(
            location = rep(this_location, 2), #ifelse(is.null(predx_point), 2, 3)),
            target = rep(this_target, 2), #ifelse(is.null(predx_point), 2, 3)),
            predx = list(predx_sample, predx_bin)
          ))
        } else {
          result <- predx::as.predx_df(list(
            location = rep(this_location, 3), #ifelse(is.null(predx_point), 2, 3)),
            target = rep(this_target, 3), #ifelse(is.null(predx_point), 2, 3)),
            predx = list(predx_sample, predx_bin, predx_point)
          ))
        }
        
        return(result)
      }
    )
  }
  
  res_csv <- merged_predx %>%
   dplyr::filter(predx_class %in% c("BinCat", "BinLwr", "Point", "PointCat")) %>%
   dplyr::mutate(
     team = "Kernel of Truth",
     mmwr_week = as.character(analysis_time_week),
     submission_date = Sys.Date(),
     unit = ifelse(target %in% c("Season onset", "Season peak week"), "week", "percent")
   ) %>%
   predx::export_flusight_csv()# %>%
#    mutate(
#      value = format(value, digits = 22, scientific = FALSE)
#    )

  csv_res_file <- file.path(file_path,
   "csv",
   paste0(
     "EW", ew,
     "-", year,
     "-", method,
     ".csv"))

  write.csv(res_csv, file = csv_res_file, row.names = FALSE)
  
  return(csv_res_file)
}


context("get_predx_forecasts_from_trajectory_samples")
library(dplyr)
library(purrr)
library(FluSight)
library(cdcfluview)
library(cdcfluutils)
library(predx)

#' Call cdcfluutils::get_predx_forecasts_from_trajectory_samples and then get
#' results in a format similar to FluSight::create_truth
predx_results_one_loc <- function(wILI, loc, year) {
  results <- bind_rows(
    cdcfluutils::get_predx_forecasts_from_trajectory_samples(
      trajectory_samples = matrix(wILI, nrow = 1),
      location = loc,
      targets = c("Season onset", "Season peak week", "Season peak percentage"),
      season = paste0(year, "/", year + 1),
      analysis_time_season_week = 50,
      predx_types = "Sample"
    ) %>%
      predx::export_csv() %>%
      mutate(
        location = ifelse(
          loc == "National",
          "US National",
          paste0("HHS Region ", substr(loc, 7, nchar(loc)))),
        forecast_week = NA
      ),
    purrr::map_dfr(
      10:43,
      function(forecast_season_week) {
        cdcfluutils::get_predx_forecasts_from_trajectory_samples(
          trajectory_samples = matrix(wILI, nrow = 1),
          location = loc,
          targets = c(paste0(1:4, " wk ahead")),
          season = paste0(year, "/", year + 1),
          analysis_time_season_week = forecast_season_week,
          predx_types = "Sample"
        ) %>%
          predx::export_csv() %>%
          mutate(
            location = ifelse(
              loc == "National",
              "US National",
              paste0("HHS Region ", substr(loc, 7, nchar(loc)))),
            forecast_week = cdcfluutils::season_week_to_year_week(
              forecast_season_week,
              weeks_in_first_season_year = cdcfluutils::get_num_MMWR_weeks_in_year(year)
            )
          )
      }
    )
  )
  
  return(results)
}


test_that("get_predx_forecasts_from_trajectory_samples agrees with FluSight::create_truth", {
  all_years <- unique(FluSight::past_baselines$year)
  flusight_truth <- purrr::map(all_years, ~ FluSight::create_truth(year = .))
  
  get_predx_forecasts_from_trajectory_samples_results <- purrr::map(all_years,
    function(year) {
      usflu <- cdcfluview::ilinet(region = "national", years = year) %>% 
        select(week, wILI = weighted_ili) %>%
        mutate(location = "US National")
      
      regionflu <- cdcfluview::ilinet(region = "HHS", years = year) %>% 
        transmute(
          location = gsub(" ", "", region),
          week = week,
          wILI = weighted_ili)
      
      w40_ind <- which(usflu$week == 40)
      us_wILI <- usflu[seq(from = w40_ind, to = nrow(usflu)), ]
      us_wILI <- us_wILI %>% pull(wILI)
      return(bind_rows(
        predx_results_one_loc(c(rep(NA, 9), us_wILI[seq_len(length(us_wILI) - 9)]), "National", year),
        purrr::map_dfr(unique(regionflu$location),
          function(loc) {
            region_wILI <- regionflu %>% filter(location == loc)
            w40_ind <- which(region_wILI$week == 40)
            region_wILI <- region_wILI[seq(from = w40_ind, to = nrow(region_wILI)), ]
            region_wILI <- region_wILI %>% pull(wILI)
            predx_results_one_loc(c(rep(NA, 9), region_wILI[seq_len(length(region_wILI) - 9)]), loc, year)
          })
      ))
    }
  )
  
  for(i in seq_along(all_years)) {
    temp <- flusight_truth[[i]] %>%
      left_join(
        get_predx_forecasts_from_trajectory_samples_results[[i]]
      ) %>%
      mutate(
        sample = ifelse(!is.na(as.numeric(sample)) & nchar(sample) < 3,
          paste0(sample, ".0"),
          sample)
      )
    
    # uncomment the following lines to get test to pass
    # not sure how big a deal these really are...
    # FluSight::create_truth sets peak to NA if after end of forecast timeframe
    #temp <- temp[!is.na(temp$bin_start_incl), ]
    
    # FluSight::create_truth only allows onset and peak to be week 41 or later; our code allows week 40
    #temp <- temp[!(temp$target == "Season onset" & temp$sample == "40.0"), ]
    #locations_peak_week_40 <- temp$location[temp$target == "Season peak week" & temp$sample == "40.0"]
    #temp <- temp[!(temp$target %in% c("Season peak week", "Season peak percentage") & temp$location %in% locations_peak_week_40), ]
    
    # FluSight::create_truth records one peak week in seasons where the same peak was reached in multiple weeks (after rounding)
    # our code keeps all such weeks
    #multi_peak_locations <- temp %>% count(location, target) %>% filter(target == "Season peak week", n > 1)
    #temp <- temp[!(temp$target == "Season peak week" & temp$location %in% multi_peak_locations$location), ]
    
    expect_identical(temp$bin_start_incl, temp$sample)
  }
})

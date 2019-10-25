## make a set of prospective prediction files for KDE/GAM model
## Nicholas Reich
## October 2017: created
## October 2018: updated
## October 2019: updated and ported to new package.

library(cdcfluutils)
library(doMC)

FIRST_YEAR_OF_CURRENT_SEASON <- 2019
this_season <- paste0(FIRST_YEAR_OF_CURRENT_SEASON, "/", FIRST_YEAR_OF_CURRENT_SEASON+1)

data(flu_data)

n_sims <- 100000

## for 2019-2020 season, 
## first predictions due on 10/28/2019 (EW44 == SW14)  MMWRweek("2019-10-28")
## using data posted on 10/25/2019 that includes up to EW42 == SW12
## last predictions due on 5/11/2020 (EW20 == SW 42) MMWRweek("2020-05-11")
## last predictions use data through EW18 == SW40
## first_analysis_time_season_week could be set to 15, but padding at front end

season_weeks <- 10:43
region_strings <- c("National", paste("Region", 1:10))
fit_path <- "inst/estimation/region-kde/fits/"

registerDoMC(3)

## fit 2019/2020 models
foreach(reg=region_strings) %dopar% {
    
  ## fit models on training seasons, using only prospective data, not LOSO
  ## this function call saves a set of .rds files that contain the list defining a "KDE fit" 
  ## one fit for each (prospective season, region) pair
  
  ## temp fix for zeroes in data
  zero_idx <- which(flu_data$weighted_ili==0)
  flu_data[zero_idx,"weighted_ili"] <- NA 
  
  # reg = region_strings[1]
  fit_region_kdes(flu_data, 
    region=reg,
    first_fit_year = FIRST_YEAR_OF_CURRENT_SEASON,
    last_fit_year = FIRST_YEAR_OF_CURRENT_SEASON,
    first_fit_week = 20, 
    path = fit_path)
}
## will receive the following warning when a region is fit and there exists a season with no onset:

## Warning message:
## In cbind(unique(data$season), as.numeric(unlist(observed_seasonal_quantities_by_season[,  :
##     NAs introduced by coercion

## make entry files
foreach(season_week = season_weeks) %dopar% {
  ## season_week <- season_weeks[1]
  make_one_kde_prediction_file(save_path = "inst/prospective-predictions/region-kde/",
    fits_path = fit_path,
    season = this_season,
    season_week = season_week,
    n_sim = n_sims)
}


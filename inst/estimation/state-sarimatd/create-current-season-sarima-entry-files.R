## make a set of prospective prediction files for KDE/GAM model
## Graham Gibson
## October 2017: created
## October 2018: updated
## October 2019: updated and ported to new package.

library(cdcfluforecasts)
library(sarimaTD)
library(jsonlite)
library(doMC)
library(predx)
library(cdcfluutils)
FIRST_YEAR_OF_CURRENT_SEASON <- 2019
this_season <- paste0(FIRST_YEAR_OF_CURRENT_SEASON, "/", FIRST_YEAR_OF_CURRENT_SEASON+1)

flu_data <- download_and_preprocess_state_flu_data()

n_sims <- 1000

## for 2019-2020 season, 
## first predictions due on 10/28/2019 (EW44 == SW14)  MMWRweek("2019-10-28")
## using data posted on 10/25/2019 that includes up to EW42 == SW12
## last predictions due on 5/11/2020 (EW20 == SW 42) MMWRweek("2020-05-11")
## last predictions use data through EW18 == SW40
## first_analysis_time_season_week could be set to 15, but padding at front end

season_weeks <- 10
region_strings <- unique(flu_data$region)
fit_path <- "inst/estimation/region-kcde/fits/"

registerDoMC(3)

rgn_lower_case <- c(
  'al', 'ak', 'az', 'ar', 'ca', 'co', 'ct', 'de', 'fl', 'ga', 'hi', 'id', 'il',
  'in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt',
  'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or',
  'pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy',
  'as','mp', 'dc', 'gu', 'pr', 'vi', 'ord', 'lax', 'jfk')

regions <- c("Alabama","Alaska","Arizona", "Arkansas","California","Colorado","Connecticut",
             "Delaware","Florida","Georgia","Hawaii","Idaho",
             "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
             "Maine","Maryland","Massachusetts","Michigan",
             "Minnesota","Mississippi","Missouri","Montana",
             "Nebraska","Nevada","New Hampshire","New Jersey",
             "New Mexico","New York","North Carolina","North Dakota",
             "Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island",
             "South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont",
             "Virginia","Washington","West Virginia","Wisconsin","Wyoming","As","mp","District of Columbia",
             "gu","Puerto Rico","Virgin Islands","ord","lax","New York City")

predx_list <- list()
idx <- 1
## fit 2019/2020 models
#foreach(reg=region_strings) %dopar% {
for (reg in region_strings){ 
  if (reg != "Commonwealth of the Northern Mariana Islands" & reg!= "Florida" ){
    ## fit models on training seasons, using only prospective data, not LOSO
    ## this function call saves a set of .rds files that contain the list defining a "KDE fit" 
    ## one fit for each (prospective season, region) pair
    
    reg_lower_case_idx <- which(regions==reg)
    cur_reg_lower_case <- rgn_lower_case[reg_lower_case_idx]
    
    cur_reg_upper_case <- regions[reg_lower_case_idx]
    
    
   # saveRDS(kcde_fit,paste0("./inst/estimation/state-kcdetd/fits/kcde_fit_",cur_reg_lower_case))
    
    analysis_time_season = "2019/2020"
    analysis_time_season_week <- flu_data$season_week[nrow(flu_data)]
    weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(analysis_time_season)
    analysis_time_ind <- nrow(flu_data)
    last_analysis_time_season_week = 41
    max_prediction_horizon <- max(4L, last_analysis_time_season_week + 
                                    1 - analysis_time_season_week)
    first_season_obs_ind <- min(which(flu_data$season == analysis_time_season))
    first_analysis_time_season_week = 10
    
    sarimaFit <- sarimaTD::fit_sarima(tail(flu_data[flu_data$region == cur_reg_upper_case,]$unweighted_ili,300),
                         ts_frequency = 52)
    
    preds <-    simulate(
        object = sarimaFit,
        nsim = 1000,
        seed = 1,
        newdata = flu_data[flu_data$region == cur_reg_upper_case,]$unweighted_ili,
        h = max_prediction_horizon
      )
    
    if (nchar(tail(flu_data$week,1)) == 1){
      current_season_epiweek <- paste0(substr(tail(flu_data$season,1),6,10),"0",tail(flu_data$week,1))
    } else{
      current_season_epiweek <- paste0(substr(tail(flu_data$season,1),1,4),tail(flu_data$week,1))
    }
    
    req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",cur_reg_lower_case,"&epiweeks=",KCDETD::move_k_week_ahead(current_season_epiweek,1)))
    nowcast_json <- jsonlite::prettify(rawToChar(req$content))
    nowcast_obj_1_wk_ahead <- fromJSON(nowcast_json)
    oneweek_ahead_nowcast <- nowcast_obj_1_wk_ahead$epidata$value
    
    req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",cur_reg_lower_case,"&epiweeks=",KCDETD::move_k_week_ahead(current_season_epiweek,2)))
    nowcast_json <- jsonlite::prettify(rawToChar(req$content))
    nowcast_obj_2_wk_ahead <- fromJSON(nowcast_json)
    twoweek_ahead_nowcast <- nowcast_obj_2_wk_ahead$epidata$value
    
    if (!is.null(oneweek_ahead_nowcast)){
      preds[,1] <- .5*preds[,1] + .5*oneweek_ahead_nowcast
    }
    
    season_start_epiweek <- 30
    if (tail(flu_data$week,1) <= 20){
      time_in <- cdcfluutils::get_num_MMWR_weeks_in_first_season_year(season) - season_start_epiweek + tail(flu_data$week,1)
    } else{
      time_in <- tail(flu_data$week,1) - season_start_epiweek + 1
    }
    
    
    sampled_historical <-rRevisedILI_fast(n = 1000,tail(flu_data[flu_data$region==cur_reg_upper_case,]$unweighted_ili,time_in),region = cur_reg_lower_case,add_nowcast = FALSE,epiweek_idx = tail(flu_data$week,1),season=tail(flu_data$season,1))
    
    preds <- cbind(preds,sampled_historical)
    preds[preds <0 ] <- 0
    predx_list[[idx]] <- get_predx_forecasts_from_trajectory_samples(trajectory_samples = preds, 
                                                                     location = cur_reg_upper_case, targets = c("Season peak week", 
                                                                                                                "Season peak percentage", paste0(1:4, " wk ahead")), 
                                                                     season = analysis_time_season, analysis_time_season_week = analysis_time_season_week, 
                                                                     first_analysis_time_season_week = first_analysis_time_season_week, 
                                                                     last_analysis_time_season_week = last_analysis_time_season_week, 
                                                                     predx_types = c("Sample", "Bin", "Point"))
    idx <- idx + 1
  }
}
library(data.table)
pred_to_write <- rbind_list(predx_list) 

submission_df <- predx_to_submission_df(pred_to_write, ew = tail(flu_data$week,1), year = substr(tail(flu_data$season,1),1,4), team = "Kernel of Truth")
## will receive the following warning when a region is fit and there exists a season with no onset:
write.csv(submission_df, file =paste0("inst/prospective-predictions/state-sarimatd/EW",tail(flu_data$week,1),"-",substr(tail(flu_data$season,1),6,10),"-ReichLab_sarimatd.csv"))
FluSight::verify_entry(submission_df,challenge ="state_ili" )
## Warning message:
## In cbind(unique(data$season), as.numeric(unlist(observed_seasonal_quantities_by_season[,  :
##     NAs introduced by coercion




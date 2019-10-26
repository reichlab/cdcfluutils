## code for KDE prediction models
## Nicholas Reich
## 7 October 2016
## 15 June 2017 - updated for cdcFlu20172018 package 

## General strategy

## ESTIMATION
## write a function 
##   inputs: a dataset
##   outputs: a list of KDE fits, one for each target
## iterate function across all regions and years to create a set of fits 
##    to all data to up the current year (not LOSO)
## 
## notes by target:
##     weekly incidence: use observations from target week +/- 1 week 
##     peak incidence: truncate above previous observations?
##     onset and peak week: discrete, but use continuous KDEs and round
##     onset week: there is a none category, consider including a point-mass based on past seasons not above baseline
##     
## PREDICTION
## write a predict function 
##   input: a fitted object from above 
##   input: a time at which to make predictions for (data not needed!)
##   input: n_sims
##   output: data_frame in format as others, specifying predictive distribution
## iterate across regions, left-out years to create predictions
## 


##' Fits and saves KDE models for all region-year combinations
##' 
##' @param data cleaned usflu dataset
##' @param region character string identifying a region
##' @param first_fit_year along with first_fit_week, defines the first time for which a fit is made
##' @param last_fit_year last year in which a fitted season begins
##' @param first_fit_week see first_fit_week, week defined on calendar scale
##' @param path filepath for saving 
##' 
##' @return nothing, just saving files
##' 
##' @details This function call assumes that the data object contains no data 
##' from the testing phase. Starting with first_fit_year/first_fit_week, a 
##' prospective fit will be made for each season until the end of the dataset.
##' 
##' @export
##' 
fit_region_kdes <- function(data, region, first_fit_year, last_fit_year, first_fit_week, path) {
    require(MMWRweek)
    require(dplyr)
    
    ### get proportion of region-seasons (across all regions and all seasons
    ### before first fit year/week) with no onset
    onsets_by_region_season <- NULL
    for(region_val in unique(data$region)) {
        
        data_subset <- data[which(data$region == region_val),]

        observed_seasonal_quantities_by_season <- t(sapply(
            as.character(unique(data_subset$season)),
            function(season) {
                get_observed_seasonal_quantities(
                    data = data_subset,
                    season = season,
                    first_CDC_season_week = 10,
                    last_CDC_season_week = 41,
                    onset_baseline = 
                        get_onset_baseline(region = region_val, season = season),
                    incidence_var = "weighted_ili",
                    incidence_bins = data.frame(
                        lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
                        upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
                    incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
                )
            }
        ))
        
        onsets_by_region_season <- bind_rows(onsets_by_region_season,
            data.frame(
                region = region_val,
                season = rownames(observed_seasonal_quantities_by_season),
                onset = as.character(unname(unlist(observed_seasonal_quantities_by_season[, 1]))),
                stringsAsFactors = FALSE
            )
        )
    }
    onsets_by_region_season$season_year <- as.numeric(substr(onsets_by_region_season$season, 1, 4))
    
    ### subsetting data

    ## subset to region of interest
    dat <- data[which(data$region == region),]
    
    ## check to make sure that all dates are in order!!
    if(any(as.numeric(diff(dat$time))<0))
        stop("dates in flu data are not ordered within region")
    
    ## assumes region is either "National" or "Region k" format
    reg_string <- ifelse(region=="National", "National", gsub(" ", "", region))
    
    ### loop over and fit each season of data, starting with first_fit_year
    years_to_fit <- first_fit_year:last_fit_year
    
    for(year_to_fit in years_to_fit) {
        ## subset to include only times before current year_to_fit
        idx <- which(dat$year == year_to_fit & dat$week == first_fit_week)
        tmpdat <- dat[seq_len(idx - 1), , drop = FALSE]

        ### create filename for saving and check to see if fit exists already
        season_to_fit <- paste0(year_to_fit, "/", (year_to_fit+1))
        filename <- paste0(
            path,
            "kde-",
            reg_string,
            "-fit-prospective-",
            gsub("/", "-", season_to_fit),
            ".rds")
        if(file.exists(filename)){
            message(paste(region, season_to_fit, "already fit ::", Sys.time()))
            next
        }
        
        ### create fits
        kde_onset_week <- fit_kde_onset_week(tmpdat,
          prob_no_onset = mean(onsets_by_region_season$onset[onsets_by_region_season$season_year < year_to_fit] == "none"))
        kde_peak_week <- fit_kde_peak_week(tmpdat)
        kde_log_peak_week_inc <- fit_kde_log_peak_week_inc(tmpdat)
        kde_log_weekly_inc <- fit_kde_weekly_inc(tmpdat)
            
        kde_fits <- list(onset_week = kde_onset_week, 
                         peak_week = kde_peak_week, 
                         log_peak_week_inc = kde_log_peak_week_inc, 
                         log_weekly_inc = kde_log_weekly_inc)
        
        ### save fits
        saveRDS(kde_fits, file = filename)
    }
}


#' Estimate KDE for onset week
#'
#' @param data data with one season left out
#' @param prob_no_onset probability to assign to "no onset" bin
#' 
#' @export
#'
#' @return a fit from density(), along with points to evaluate at
fit_kde_onset_week <- function(data, prob_no_onset) {
    require(dplyr)
    ### get onset
    observed_seasonal_quantities_by_season <- t(sapply(
        unique(data$season),
        function(season) {
            get_observed_seasonal_quantities(
                data = data,
                season = season,
                first_CDC_season_week = 10,
                last_CDC_season_week = 41,
                onset_baseline = 
                    get_onset_baseline(region = data$region[1], season = season),
                incidence_var = "weighted_ili",
                incidence_bins = data.frame(
                    lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
                    upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
                incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
            )
        }
    ))
    
    ## vector of onset weeks
    onset_week <- as.numeric(unlist(observed_seasonal_quantities_by_season[, 1]))
    
    ### calculate density based on non-NA values
    onset_week <- onset_week[!is.na(onset_week)]
    return(list(kde=density(onset_week, na.rm=TRUE), prob_no_onset=prob_no_onset, x=onset_week))
}


#' Estimate KDE for peak week
#'
#' @param data 
#'
#' @return a fit from density(), along with points to evaluate at 
fit_kde_peak_week <- function(data) {
    require(dplyr)
    ### get peaks
    peak_week <- data %>% group_by(season) %>%
        ## filter to ensure estimated peaks are within CDC-specified forecasting season
        filter(season_week >= 10 & season_week <= 42) %>%
        summarize(peak_week = season_week[which(weighted_ili == max(weighted_ili, na.rm=TRUE))]) %>%
        .$peak_week

    ### calculate density
    return(list(kde=density(peak_week, na.rm=TRUE), x=peak_week))
}

#' Estimate KDE for LOG peak week incidence
#'
#' @param data 
#'
#' @return a fit from density(), along with points to evaluate at
fit_kde_log_peak_week_inc <- function(data) {
    require(dplyr)
    ### get peaks
    peak_week_inc <- data %>% group_by(season) %>%
        summarize(peak_week_inc = weighted_ili[which(weighted_ili == max(weighted_ili, na.rm=TRUE))]) %>%
        .$peak_week_inc
    
    ### calculate density
    return(list(kde=density(log(peak_week_inc), na.rm=TRUE), x=log(peak_week_inc)))
}

#' Fit a GAM for weekly incidence
#'
#' @param data 
#' 
#' @export
#'
#' @return a mgcv smooth spline fit
fit_kde_weekly_inc <- function(data) {
    require(mgcv)
    
    ### fit model
    fm <- gam(log(weighted_ili) ~ s(season_week, bs = "cc"),
              data=data)
    
    return(fm)
}



#' Wrapper for prediction and evaluation of KDE fits
#'
#' @param data data, to be used for evaluation
#' @param region region of focus
#' @param path filepath to fitted models
#' @param n_sim number of simulations to run for predictive distributions
#'
#' @return NULL just saves a file
#'
#' @description The function follows this outline
#' \itemize{
#'       \item{determine set of years for which LOSO fits are available}
#'       \item{For those years, generate predictions made at each season-week. 
#'             These will be the same predictions repeated, as models not dependent on any inputs.}
#'       \item{Compare predictions to the real data}
#'       \item{Save object with comparisons}
#' }
#' 
#' @export
predict_region_kde <- function(data, region, path, n_sim) {
    
    require(dplyr)
    
    ### SETUP
    ## load baselines for onset calculations
    baselines <- read.csv(file="data-raw/cdc-baselines.csv")
    ## NOTE: assumes region is either "X" or "Region k" format
    reg_string <- ifelse(region=="X", "National", gsub(" ", "", region))
    idx <- which(baselines$region==reg_string & baselines$season=="2015/2016")
    reg_baseline <- baselines[idx, "baseline"]
    
    ## subset data to be only the region of interest
    dat <- data[which(data$region == region),]
    
    ### SET GLOBAL VALUES

    ## incidence bins for predictions of incidence in individual weeks and at the peak
    inc_bins <- c(0, seq(from = .05, to = 12.95, by = 0.1), Inf)
    incidence_bin_names <- as.character(seq(from = 0, to = 13, by = 0.1))
    
    ## for 2016-2017 season, first predictions due on 11/7/2016 (EW45 == SW15)
    ## using data posted on 11/4/2016 that includes up to EW43 == SW13
    ## last predictions due on 5/15/2017 (EW20 == SW 42)
    ## last predictions use data through EW40 == SW18
    ## first_analysis_time_season_week could be set to 13, but padding at front end
    first_analysis_time_season_week <- 10 # == week 40 of year
    last_analysis_time_season_week <- 41 # == week 19 or 18, depending on whether a 53 week season
    
    ## subset of season weeks for CDC competition
    cdc_season_weeks <- as.character(10:42)
    
    ### determine which years we have fits for
    fnames <- system(paste0('ls ', path), intern=TRUE)
    fnames_this_reg <- fnames[grep(paste0(reg_string,"-"), fnames)] ## need hyphen for Reg 1 vs. 10
    where_is_the_dot <- regexpr(pattern = "\\.", fnames_this_reg)
    seasons_this_reg <- substr(fnames_this_reg, start=where_is_the_dot - 9, stop=where_is_the_dot-1)
    seasons_this_reg <- gsub("-", "/", seasons_this_reg)
    
    ### calculate observed values of targets across seasons for later comparison
    ## season week is week of flu season (week 1 is in Oct)
    ## week is epidemic week of calendar year (week 1 is in January)
    
    ## onset week
    onset_weeks <- dat %>% group_by(season) %>% 
        mutate(ili_lag1 = lag(weighted_ili, 1),
               ili_lag2 = lag(weighted_ili, 2),
               ili_lag3 = lag(weighted_ili, 3),
               onset = ili_lag1 > reg_baseline & ili_lag2 > reg_baseline & ili_lag3 > reg_baseline) %>%
        filter(onset) %>%
        summarise(season_week = min(season_week)-2) %>% ## may not do the right thing if onset is in first 2 weeks
        left_join(dat) %>%
        select(season, season_week, week)
    
    ## peak week
    peak_weeks <- dat %>% group_by(season) %>%
        summarize(season_week = season_week[which(weighted_ili == max(weighted_ili, na.rm=TRUE))]) %>%
        left_join(dat) %>%
        select(season, season_week, week)
    
    ## peak week incidence
    peak_week_inc <- dat %>% group_by(season) %>%
        summarize(peak_week_inc = weighted_ili[which(weighted_ili == max(weighted_ili, na.rm=TRUE))])
    
    ### generate predictions made at each season-week
    ## make dataset for storage
    predictions_df <- make_predictions_dataframe(dat, 
                                                 model_name="kde", 
                                                 incidence_bin_names=incidence_bin_names)
    
    results_save_row <- 1L
    
    for(analysis_time_season in seasons_this_reg) {
        message(paste(region, analysis_time_season, "::", Sys.time()))
        ## Only do something if there is something to predict in the season that would be held out
        if(!all(is.na(dat$weighted_ili[dat$season == analysis_time_season]))) {
            ## load KDE fit
            kde_fit <- readRDS(file = paste0(
                path,
                "kde-",
                reg_string,
                "-fit-leave-out-",
                gsub("/", "-", analysis_time_season),
                ".rds"))
            
            ## make predictions for each prediction target in the left-out season
            ## for each possible "last observed" week, starting with the last week of the previous season
            first_season_ind <- min(which(dat$season == analysis_time_season))
            last_season_ind <- max(which(dat$season == analysis_time_season))
            
            last_analysis_time_season_week_in_dat <- max(dat$season_week[dat$season == analysis_time_season])
            
            ## make prediction for this analysis_time_season
            onset_week_preds <- predict_kde_onset_week(kde_fit$onset_week, n_sim)
            peak_week_preds <- predict_kde_peak_week(kde_fit$peak_week, n_sim)
            peak_week_inc_preds <- predict_kde_log_peak_week_inc(kde_fit$log_peak_week_inc, 
                                                                 bins=inc_bins,
                                                                 bin_names=incidence_bin_names, 
                                                                 n_sim)
            weekly_inc_preds <- predict_kde_log_weekly_inc(fm = kde_fit$log_weekly_inc, 
                                                           season_weeks = 1:53, 
                                                           bins = inc_bins, 
                                                           bin_names = incidence_bin_names,
                                                           n_sim = n_sim)
            
            ## calculate observed values to compare against
            observed_onset_week <- unlist(onset_weeks[which(onset_weeks$season==analysis_time_season),"season_week"])
            observed_peak_week <- unlist(peak_weeks[which(peak_weeks$season==analysis_time_season),"season_week"])
            observed_peak_week_inc <- unlist(peak_week_inc[which(peak_week_inc$season==analysis_time_season),"peak_week_inc"])
            
            ### calculate log score for unchanging predictions
            
            ## calculate onset week log scores
            if(length(observed_onset_week)==0) { ## if no observed onset, assign none
                obs_onset_week_char <- "none"
                log_score_for_onset_week <- compute_competition_log_score(log(onset_week_preds[c(cdc_season_weeks, "none")]),
                                                                          observed_bin = obs_onset_week_char,
                                                                          prediction_target = "onset_week")
            } else if(observed_onset_week<10) { ## if onset prior to week 10, assign log-score of NA
                log_score_for_onset_week <- NA
            } else if(observed_onset_week>42) { ## if onset after week 42, assign "none"
                obs_onset_week_char <- "none"
                log_score_for_onset_week <- compute_competition_log_score(log(onset_week_preds[c(cdc_season_weeks, "none")]),
                                                                          observed_bin = obs_onset_week_char,
                                                                          prediction_target = "onset_week")
            } else { ## otherwise, onset week in [10, 42]
                obs_onset_week_char <- as.character(observed_onset_week)
                log_score_for_onset_week <- compute_competition_log_score(log(onset_week_preds[c(cdc_season_weeks, "none")]),
                                                                          observed_bin = obs_onset_week_char,
                                                                          prediction_target = "onset_week")
            }

            ## calculate peak week log scores
            ## assign the index based on output from predict_kde_peak 
            log_score_for_peak_week <- ifelse(observed_peak_week<10 | observed_peak_week > 42,
                                              NA,
                                              compute_competition_log_score(log(peak_week_preds[cdc_season_weeks]),
                                                                     observed_bin = as.character(observed_peak_week),
                                                                     prediction_target = "peak_week"))
            
            ## calculate peak week incidence log scores
            peak_inc_bin_char <- get_inc_bin(observed_peak_week_inc)
            log_score_for_peak_week_inc <- compute_competition_log_score(log(peak_week_inc_preds),
                                                                         observed_bin = peak_inc_bin_char,
                                                                         prediction_target="peak_inc")

            ## sequence to loop over
            time_seq <- seq(from = first_analysis_time_season_week, to = min(last_analysis_time_season_week, last_analysis_time_season_week_in_dat) - 1)
            
            for(analysis_time_season_week in time_seq) {
                ## assign log scores
                predictions_df[results_save_row, "onset_log_score"] <- log_score_for_onset_week
                predictions_df[results_save_row, "peak_week_log_score"] <- log_score_for_peak_week
                predictions_df[results_save_row, "peak_inc_log_score"] <- log_score_for_peak_week_inc
                
                ## calculate and assign peak week incidence log scores
                analysis_time_ind <- which(dat$season == analysis_time_season &
                                               dat$season_week == analysis_time_season_week)

                ## 1 week ahead
                observed_ph_1_inc <- dat[analysis_time_ind+1, "weighted_ili"]
                ph_1_inc_bin <- get_inc_bin(observed_ph_1_inc)
                ph_1_season_week <- analysis_time_season_week+1
                log_score_for_ph_1_inc <- ifelse(is.na(ph_1_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_1_season_week]),
                                                                        observed_bin = ph_1_inc_bin,
                                                                        prediction_target="ph1_inc"))
                predictions_df[results_save_row, "ph_1_inc_log_score"] <- log_score_for_ph_1_inc
                
                ## 2 week ahead
                observed_ph_2_inc <- dat[analysis_time_ind+2, "weighted_ili"]
                ph_2_inc_bin <- get_inc_bin(observed_ph_2_inc)
                ph_2_season_week <- analysis_time_season_week+2
                log_score_for_ph_2_inc <- ifelse(is.na(ph_2_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_2_season_week]),
                                                                        observed_bin = ph_2_inc_bin,
                                                                        prediction_target="ph2_inc"))
                predictions_df[results_save_row, "ph_2_inc_log_score"] <- log_score_for_ph_2_inc
                
                ## 3 week ahead
                observed_ph_3_inc <- dat[analysis_time_ind+3, "weighted_ili"]
                ph_3_inc_bin <- get_inc_bin(observed_ph_3_inc)
                ph_3_season_week <- analysis_time_season_week+3
                log_score_for_ph_3_inc <- ifelse(is.na(ph_3_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_3_season_week]),
                                                                        observed_bin = ph_3_inc_bin,
                                                                        prediction_target="ph3_inc"))
                predictions_df[results_save_row, "ph_3_inc_log_score"] <- log_score_for_ph_3_inc
                
                ## 4 week ahead
                observed_ph_4_inc <- dat[analysis_time_ind+4, "weighted_ili"]
                ph_4_inc_bin <- get_inc_bin(observed_ph_4_inc)
                ph_4_season_week <- analysis_time_season_week+4
                log_score_for_ph_4_inc <- ifelse(is.na(ph_4_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_4_season_week]),
                                                                        observed_bin = ph_4_inc_bin,
                                                                        prediction_target="ph4_inc"))
                predictions_df[results_save_row, "ph_4_inc_log_score"] <- log_score_for_ph_4_inc
                
                ## set other fixed stuff
                predictions_df[results_save_row, "analysis_time_season"] <- analysis_time_season
                predictions_df[results_save_row, "analysis_time_season_week"] <- analysis_time_season_week
                predictions_df[results_save_row, "prediction_week_ph_1"] <- analysis_time_season_week + 1
                predictions_df[results_save_row, "prediction_week_ph_2"] <- analysis_time_season_week + 2
                predictions_df[results_save_row, "prediction_week_ph_3"] <- analysis_time_season_week + 3
                predictions_df[results_save_row, "prediction_week_ph_4"] <- analysis_time_season_week + 4
                
                results_save_row <- results_save_row + 1
                
            } # analysis_time_season_week
        }
    } # analysis_time_season
    
    ## if there are extra rows in the predictions_df, delete them
    if(results_save_row <= nrow(predictions_df)) {
        predictions_df <- predictions_df[
            -seq(from = results_save_row, to = nrow(predictions_df)),
            ,
            drop = FALSE
            ]
    }
    
    ### save and return
    kde_predictions_df <- predictions_df
    filename <- paste0(
        "inst/estimation/loso-predictions/kde-",
        reg_string,
        "-loso-predictions.rds")
    saveRDS(kde_predictions_df, file = filename)
    
}

#' Wrapper for prediction of a single time-point from KDE fits
#'
#' @param fits_path filepath to fitted models
#' @param save_path filepath to save models in
#' @param season character string of season for current fit, in format "XXXX/YYYY"
#' @param season_week numeric season week of first week for which predictions are desired
#' @param n_sim number of simulations to run for predictive distributions
#'
#' @return NULL just saves a file
#' @export
#'
make_one_kde_prediction_file <- function(fits_path, save_path, season, season_week, n_sim){
    require(MMWRweek)
    ## determine if a 53-week season
    last_day_cal_year <- as.Date(paste0(substr(season, 0, 4),"-12-31"))
    season_has_EW53 <- MMWRweek(last_day_cal_year)$MMWRweek==53
    num_EW <- ifelse(season_has_EW53, 53, 52)
    
    ## set globals for incidence bins and names
    inc_bins <- c(0, seq(from = .05, to = 12.95, by = 0.1), Inf)
    incidence_bin_names <- as.character(seq(from = 0, to = 13, by = 0.1))
    
    region_strings <- c("National", paste0("Region", 1:10))
    cdc_region_strings <- c("US National", paste("HHS Region", 1:10))
    
    epiweek <- (season_week-1+30)%%num_EW + 1
    epiweek_year <- ifelse(epiweek<35, ## if epiweek is low, we are in end of season
        substr(season, 6, 9), 
        substr(season, 1, 4))
    fname <- paste0(save_path,"/EW", sprintf("%02d", epiweek), "-", epiweek_year, "-ReichLab_kde.csv")
    
    for(i in 1:length(region_strings)) {
        region <- region_strings[i]
        cdc_region <- cdc_region_strings[i]
        ## load KDE fit
        kde_fit <- readRDS(file = paste0(
            fits_path,
            "kde-",
            region,
            "-fit-prospective-",
            gsub("/", "-", season),
            ".rds"))
        tmp <- make_kde_region_prediction(kde_fit, 
            cdc_region_string=cdc_region, 
            inc_bins=inc_bins, 
            inc_bin_names=incidence_bin_names, 
            season_has_EW53=season_has_EW53,
            season_week = season_week,
            n_sim=n_sim)
        add_colnames <- ifelse(i==1, TRUE, FALSE) ## only include rownames in first file
        write.table(tmp, file=fname, quote=FALSE, sep=",", row.names = FALSE,
            append=!add_colnames, col.names=add_colnames)
    }
}

#' Make a KDE prediction for one region
#'
#' @param fit a kde_fit R object
#' @param cdc_region_string string for one CDC region
#' @param inc_bins list of incidence bins
#' @param inc_bin_names names of incidence bins
#' @param season_has_EW53 indicator of whether the current season has a week 53
#' @param season_week string of the season week for which to make the prediction for
#' @param n_sim number of simulations to run
#'
#' @return data.frame of the prediction, in CDC format
#' @export
make_kde_region_prediction <- function(fit, cdc_region_string, inc_bins, inc_bin_names, season_has_EW53, season_week, n_sim) {
    num_EW <- ifelse(season_has_EW53, 53, 52)
    
    ## read in output template
    if(season_has_EW53) {
        out <- read.csv("data-raw/region-prediction-template-EW53.csv")
    } else {
        out <- read.csv("data-raw/region-prediction-template.csv")
    }
    
    out$Location <- cdc_region_string
    
    ### ONSETS
    onset_week_preds <- predict_kde_onset_week(fit$onset_week, n_sim)
    ## truncating to weeks observed and adding delta to ensure no -Inf log scores
    last_week_bin <- ifelse(season_has_EW53, 43, 42)
    onset_week_preds_used <- onset_week_preds[c(10:last_week_bin, 53)] + 1/n_sim
    onset_week_preds_std <- onset_week_preds_used/sum(onset_week_preds_used)
    idx_onset_bins <- which(out$Target=="Season onset" & out$Type=="Bin")
    out[idx_onset_bins, "Value"] <- onset_week_preds_std
    idx_onset_point <- which(out$Target=="Season onset" & out$Type=="Point")
    onset_season_week <- calc_median_from_binned_probs(onset_week_preds_std)
    out[idx_onset_point, "Value"] <- (onset_season_week-1+30)%%num_EW +1
    
    ### PEAKS
    peak_week_preds <- predict_kde_peak_week(fit$peak_week, n_sim)
    ## truncating to weeks observed and adding delta to ensure no -Inf log scores
    peak_week_preds_used <- peak_week_preds[10:last_week_bin] + 1/n_sim
    peak_week_preds_std <- peak_week_preds_used/sum(peak_week_preds_used)
    idx_peak_bins <- which(out$Target=="Season peak week" & out$Type=="Bin")
    out[idx_peak_bins, "Value"] <- peak_week_preds_std
    idx_peak_point <- which(out$Target=="Season peak week" & out$Type=="Point")
    peak_season_week <- calc_median_from_binned_probs(peak_week_preds_std)
    out[idx_peak_point, "Value"] <- (peak_season_week-1+30)%%num_EW +1
    
    ### PEAK INC
    peak_inc_preds <- predict_kde_log_peak_week_inc(fit$log_peak_week_inc, 
        bins=inc_bins,
        bin_names=inc_bin_names, 
        n_sim)
    peak_inc_preds_used <- peak_inc_preds + 1/n_sim
    peak_inc_preds_std <- peak_inc_preds_used/sum(peak_inc_preds_used)
    idx_peak_inc_bins <- which(out$Target=="Season peak percentage" & out$Type=="Bin")
    ## check that names match
    if(any(names(peak_inc_preds) != out[idx_peak_inc_bins, "Bin_start_incl"])){
        error("Peak incidencebin names don't match.")
    }
    out[idx_peak_inc_bins, "Value"] <- peak_inc_preds_std
    idx_peak_inc_point <- which(out$Target=="Season peak percentage" & out$Type=="Point")
    out[idx_peak_inc_point, "Value"] <- calc_median_from_binned_probs(peak_inc_preds_std)
    
    
    ## WEEKLY INCIDENCE
    weekly_inc_preds <- predict_kde_log_weekly_inc(fm = fit$log_weekly_inc, 
        season_weeks = season_week:(season_week+3), 
        bins = inc_bins, 
        bin_names = inc_bin_names,
        n_sim = n_sim)
    weekly_inc_preds <- weekly_inc_preds + 1/n_sim
    weekly_inc_preds <- apply(weekly_inc_preds, 
        FUN=function(x) x/sum(x),
        MAR=2)
    for(i in 1:4){
        idx_weekly_bins <- which(out$Target==paste(i, "wk ahead") & out$Type=="Bin")
        out[idx_weekly_bins, "Value"] <- weekly_inc_preds[,i]
        idx_weekly_inc_point <- which(out$Target==paste(i, "wk ahead") & out$Type=="Point")
        weekly_inc_preds_vec <- weekly_inc_preds[,i]
        names(weekly_inc_preds_vec) <- inc_bin_names
        out[idx_weekly_inc_point, "Value"] <- calc_median_from_binned_probs(weekly_inc_preds_vec)
    }
    
    return(out)
}

#' Functions for predicting from KDE fits
#'


#' Compute predictive distribution for onset_week from KDE fit
#'
#' @param fm output from fit_kde_onset_week()
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return vector of 53 probabilities representing the predictive distribution of onset week, 
#'         in bins representing season week 1 through 52 (or EW31 through EW30 -- or EW29 if a 53 week year) 
#'         with the final entry representing the probability of no onset
#'
#' @export
predict_kde_onset_week <- function(fm, n_sim) {
    kde_onset_week <- fm$kde
    prob_no_onset <- fm$prob_no_onset
    onsets <- fm$x
    x.new <- rnorm(n_sim, 
        mean=sample(onsets, size = n_sim, replace = TRUE), 
        sd=kde_onset_week$bw)
    pred_onsets_rounded <- round(x.new)
    
    ## removing early/late onset draws
    idx_out_of_bounds <- which(pred_onsets_rounded<=0 | pred_onsets_rounded>52)
    if(length(idx_out_of_bounds)>0){
        pred_onsets_rounded <- pred_onsets_rounded[-idx_out_of_bounds] 
    }
    
    onsets_binned <- tabulate(pred_onsets_rounded, nbins=52)
    ## calculate number of "nones" to pad for correct probability
    nones_to_add <- round(prob_no_onset/(1-prob_no_onset)*length(pred_onsets_rounded))
    
    onset_bin_probs <- c(onsets_binned, nones_to_add)/(nones_to_add+length(pred_onsets_rounded))
    names(onset_bin_probs) <- c(1:52, "none")
    
    return(onset_bin_probs)    
}

#' Compute predictive distribution for peak week from KDE fit
#'
#' @param fm output from fit_kde_peak_week()
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return vector of 52 probabilities representing the predictive distribution of peak week, 
#'         in bins representing season week 1 through 52 (or EW31 through EW30 -- or EW29 if a 53 week year) 
#'         with the final entry representing the probability of no onset#'
#'
#' @export
predict_kde_peak_week <- function(fm, n_sim) {
    kde_peak_week <- fm$kde
    peaks <- fm$x
    x.new <- rnorm(n_sim, 
        mean=sample(peaks, size = n_sim, replace = TRUE), 
        sd=kde_peak_week$bw)
    pred_peaks_rounded <- round(x.new)
    
    ## removing early/late peak samples
    idx_out_of_bounds <- which(pred_peaks_rounded<=0 | pred_peaks_rounded>52)
    if(length(idx_out_of_bounds)>0){
        pred_peaks_rounded[-idx_out_of_bounds] 
    }
    
    peaks_binned <- tabulate(pred_peaks_rounded, nbins=52)
    
    peak_bin_probs <- peaks_binned/(n_sim-length(idx_out_of_bounds))
    names(peak_bin_probs) <- as.character(1:52)
    
    return(peak_bin_probs)    
}


#' Compute predictive distribution for peak week incidence from KDE fit
#'
#' @param fm output from fit_kde_log_peak_week_inc()
#' @param bins cutpoints for incidence bins 
#' @param bin_names vector of bin names, length 1 fewer than length(bins)
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return vector of probabilities representing the predictive distribution of NOT LOG peak week incidence, 
#'         in given bins 
#'
#' @export
predict_kde_log_peak_week_inc <- function(fm, bins, bin_names, n_sim) {
    kde_peak_week_inc <- fm$kde
    log_peak_inc <- fm$x
    x_new <- rnorm(n_sim, 
        mean=base::sample(log_peak_inc, size = n_sim, replace = TRUE), 
        sd=kde_peak_week_inc$bw)
    pred_peaks_binned <- cut(exp(x_new), breaks=bins, right=FALSE) ## CDC incidence targets are [a,b)
    peak_inc_bin_probs <- table(pred_peaks_binned)/n_sim
    
    names(peak_inc_bin_probs) <- bin_names
    return(peak_inc_bin_probs)    
}


#' Compute predictive distribution for weekly incidence from GAM fit
#'
#' @param fm a GAM fit for weekly incidence
#' @param season_weeks season_weeks to predict for
#' @param bins cutpoints for incidence bins 
#' @param bin_names vector of bin names, length 1 fewer than length(bins)
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return matrix of probabilities representing the predictive distribution of NOT LOG scale weekly incidence, 
#'         with row correponding to the bin and column corresponding to each given season_week
#'
#' @export
predict_kde_log_weekly_inc <- function(fm, season_weeks, bins, bin_names, n_sim) {
    require(mgcv)
    
    ## code adapted from example at: https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html
    
    ## extract parameter estiamtes and cov matrix...
    beta <- coef(fm)
    Vb <- vcov(fm)
    
    ## simulate replicate beta vectors from posterior...
    Cv <- chol(Vb)
    nb <- length(beta)
    br <- t(Cv) %*% matrix(rnorm(n_sim*nb),nb,n_sim) + beta
    
    ## turn these into replicate linear predictors...
    Xp <- predict(fm,newdata=data.frame(season_week=season_weeks),type="lpmatrix")
    lp <- Xp%*%br
    fv <- lp ## ... finally, replicate expected value vectors
    
    ## now simulate from normal deviates with mean as in fv
    ## and estimated scale...
    tmp <- matrix(rnorm(length(fv), mean=fv, sd=sqrt(fm$sig2)), nrow=nrow(fv), ncol(fv))
    
    #plot(rep(xp,n_sim),exp(tmp),pch=".", ylim=c(-1, 10)) ## plotting replicates
    #points(data$season_week,data$weighted_ili,pch=19,cex=.5) ## and original data
    
    ## compute 95% prediction interval...
    #PI <- apply(exp(tmp),1,quantile,prob=c(.025,0.975))
    
    weekly_inc_bin_probs <- apply(tmp, 
        MARGIN=1,
        FUN = function(x) {
            binned_inc <- cut(exp(x), breaks=bins, right=FALSE)
            table(binned_inc)/n_sim
        })
    rownames(weekly_inc_bin_probs) <- bin_names
    
    return(weekly_inc_bin_probs)    
}


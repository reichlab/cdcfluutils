# R code for SIR fit
library(readr)
library(dplyr)
library(ggplot2)
library(rstan)
library(doParallel)
library(hforecast)
library(here)
library(grid)
library(tidyr)
library(cdcfluutils)
library(predx)
library(gridExtra)

data <- readRDS("inst/make-submission-files/make-GPthingy-submission-file/us_flu.rds")
#download_and_preprocess_flu_data()

season_week <- data %>% filter(!is.na(weighted_ili)) %>% tail(1) %>% pull(season_week)

sample_files <- Sys.glob(paste0("inst/make-submission-files/make-GPthingy-submission-file/samples/samples_depth12_chain*_last_obs_2019_2020_", season_week, "_logit_normal.csv"))
#sample_files <- Sys.glob(paste0("inst/make-submission-files/make-GPthingy-submission-file/samples/samples_depth12_chain*_last_obs_2019_2020_", season_week, "_logit_t.csv"))

#sample_files <- Sys.glob("inst/make-submission-files/make-GPthingy-submission-file/samples/samples_depth12_chain*_last_obs_2019_2020_15_logit_normal_all_seasons.csv")
#sample_files <- Sys.glob("inst/make-submission-files/make-GPthingy-submission-file/samples/samples_depth12_chain*_last_obs_2019_2020_15_logit_t.csv")

#sample_files <- paste0(
#  "inst/make-submission-files/make-GPthingy-submission-file/samples/samples_depth12_chain",
#  1:2,
#  "_last_obs_2019_2020_52.csv")

# Seasons to be used as part of training or prediction data set
predict_season <- "2019/2020"
predict_last_obs_season_week <- data %>% filter(!is.na(weighted_ili)) %>% tail(1) %>% pull(season_week)
#seasons <- paste0(1997:(as.integer(substr(predict_season, 1, 4))), "/", 1998:(as.integer(substr(predict_season, 6, 9))))
seasons <- paste0(2002:(as.integer(substr(predict_season, 1, 4))), "/", 2003:(as.integer(substr(predict_season, 6, 9))))

sample_ind <- 1L
stan_data <- cdcfluutils::preprocess_reg_data_for_GPthingy_stan(
  us_flu = readRDS("inst/make-submission-files/make-GPthingy-submission-file/us_flu.rds"),
  seasons = seasons,
  predict_season = predict_season,
  predict_last_obs_season_week = predict_last_obs_season_week,
  get_all_region_seasons = TRUE,
  pred_all_t = FALSE
)

us_flu <- readRDS("inst/make-submission-files/make-GPthingy-submission-file/us_flu.rds")
unique_regions <- unique(us_flu$region)

su_with_pred <- which(stan_data$num_pred_ili > 0)
s_with_pred <- stan_data$SU_to_S_map[su_with_pred]
s_names_with_pred <- unique_regions[s_with_pred]

sample_base_by_s <- cumsum(stan_data$num_pred_ili[su_with_pred]) - stan_data$num_pred_ili[su_with_pred[1]]
sample_inds_by_s <- lapply(seq_along(su_with_pred),
                           function(i) {
                             su_ind <- su_with_pred[i]
                             sample_base_by_s[i] + seq_len(stan_data$num_pred_ili[su_ind])
                           })

all_samples <- purrr::map_dfr(
  seq_along(sample_files),
  function(chain_num) {
    cat(chain_num)
    cat("\n")
    
    sample_file <- sample_files[chain_num]
    #stan_model_samples <- rstan::read_stan_csv(sample_file)
    #temp <- rstan::extract(stan_model_samples)
    sample_lines <- read_csv(sample_file, skip = 38)
    ili_t_concat_pred_samples <- sample_lines %>% dplyr::select(starts_with("ili_t_concat_pred"))
    
    #samples_by_s <- purrr::map_dfr(
    #  seq_along(su_with_pred),
    #  function(s_ind) {
    #    samples_matrix <- ili_t_concat_pred_samples[, sample_inds_by_s[[s_ind]], drop = FALSE] - 0.00049
    #    return(samples_matrix)
    #samples_vec <- t(samples_matrix)
    #dim(samples_vec) <- prod(dim(samples_vec))
    #result <- data.frame(
    #  chain_num = chain_num,
    #  region = s_names_with_pred[s_ind],
    #  season_week = rep(seq_len(ncol(samples_matrix)), times = nrow(samples_matrix)),
    #  sample_num = rep(seq_len(nrow(samples_matrix)), each = ncol(samples_matrix)),
    #  ili_hat = samples_vec,
    #  stringsAsFactors = FALSE
    #)
    #return(result)
    #  }
    #)
    
    return(ili_t_concat_pred_samples)
  }
)

inds_to_drop <- apply(
  all_samples,
  1,
  function(sr) {all(is.na(sr))}
)

all_samples <- all_samples %>%
  slice(-inds_to_drop) %>%
  as.matrix()

weeks_in_first_season_year <-
  get_num_MMWR_weeks_in_first_season_year(predict_season)
last_analysis_time_season_week <- weeks_in_first_season_year - 11

res_predx <- purrr::map_dfr(
  seq_along(s_names_with_pred),
  function(region_ind) {
    region <- s_names_with_pred[region_ind]
    cat(region)
    cat(" ")
    cat(Sys.time())
    cat("\n")
    
    ## subset data to be only the region-specific data
    data <- us_flu[us_flu$region == region,]
    
    analysis_time_ind <- nrow(data)
    max_prediction_horizon <- max(4L,
                                  last_analysis_time_season_week + 1 - predict_last_obs_season_week)
    first_season_obs_ind <- min(which(data$season == predict_season))
    
    trajectory_samples <- round(
      (all_samples[, sample_inds_by_s[[region_ind]], drop = FALSE] - 0.00049)*100,
      digits = 1)
    
    forecast_predx <- get_predx_forecasts_from_trajectory_samples(
      trajectory_samples = trajectory_samples,
      location = region,
      targets = c("Season onset", "Season peak week", "Season peak percentage",
                  paste0(1:4, " wk ahead")),
      season = predict_season,
      analysis_time_season_week = predict_last_obs_season_week,
      first_analysis_time_season_week = 10,
      last_analysis_time_season_week = weeks_in_first_season_year - 11,
      predx_types = c("Sample", "Bin", "Point")
    )
    
    return(forecast_predx)
  }) %>%
  mutate(
    location = ifelse(
      location == "National",
      "US National",
      paste0("HHS ", location)
    )
  )

analysis_time_week <- us_flu %>%
  filter(!is.na(weighted_ili)) %>%
  pull(week) %>%
  tail(1)
analysis_time_year <- us_flu %>%
  filter(!is.na(weighted_ili)) %>%
  pull(year) %>%
  tail(1)

options(digits = 22, scipen = 9999)
csv_df <- cdcfluutils::predx_to_submission_df(
  res_predx,
  ew = analysis_time_week,
  year = analysis_time_year,
  team = "Kernel of Truth")

submissions_save_path <- "inst/submissions/region-GPthingy"
predx_res_file <- file.path(
  submissions_save_path,
  "predx",
  paste0(
    "EW", analysis_time_week,
    "-", analysis_time_year,
    "-KoT_GPthingy",
    ".rds"))
saveRDS(res_predx, predx_res_file)

csv_res_file <- file.path(
  submissions_save_path,
  "csv",
  paste0(
    "EW", analysis_time_week,
    "-", analysis_time_year,
    "-KoT_GPthingy",
    ".csv"))
write.csv(
  csv_df,
  csv_res_file,
  row.names = FALSE)


data <- download_and_preprocess_flu_data()

plot_save_path <- paste0(
  submissions_save_path,
  "/plots/",
  "EW", analysis_time_week,
  "-", analysis_time_year,
  "-KoT_GPthingy",
  ".pdf")
make_predictions_plots(
  preds_save_file = csv_res_file,
  plots_save_file = plot_save_path,
  data = data
)






all_samples_df <- purrr::map_dfr(
  seq_along(s_names_with_pred),
  function(region_ind) {
    region <- s_names_with_pred[region_ind]
    cat(region)
    cat(" ")
    cat(Sys.time())
    cat("\n")
    
    ## subset data to be only the region-specific data
    data <- us_flu[us_flu$region == region,]
    
    analysis_time_ind <- nrow(data)
    max_prediction_horizon <- max(4L,
                                  last_analysis_time_season_week + 1 - predict_last_obs_season_week)
    first_season_obs_ind <- min(which(data$season == predict_season))
    
    samples_vec <- t((all_samples[, sample_inds_by_s[[region_ind]], drop = FALSE] - 0.00049)*100)
    dim(samples_vec) <- prod(dim(samples_vec))
    return(
      data.frame(
        location = region,
        season_week = rep(seq_len(52), times = nrow(all_samples)),
        ili_hat = samples_vec
      )
    )
  }) %>%
  mutate(
    location = ifelse(
      location == "National",
      "US National",
      paste0("HHS ", location)
    )
  )


sample_summaries <- all_samples_df %>%
  group_by(location, season_week) %>%
  summarize(
    median_ili_hat = median(ili_hat - 0.049, na.rm = TRUE),
    #    q0.005 = quantile(ili_hat - 0.049, 0.005, na.rm = TRUE),
    #    q0.995 = quantile(ili_hat - 0.049, 0.995, na.rm = TRUE),
    q0.005 = quantile(ili_hat - 0.049, 0.01, na.rm = TRUE),
    q0.995 = quantile(ili_hat - 0.049, 0.99, na.rm = TRUE),
    q0.025 = quantile(ili_hat - 0.049, 0.025, na.rm = TRUE),
    q0.975 = quantile(ili_hat - 0.049, 0.975, na.rm = TRUE),
    q0.10 = quantile(ili_hat - 0.049, 0.10, na.rm = TRUE),
    q0.90 = quantile(ili_hat - 0.049, 0.90, na.rm = TRUE),
    q0.25 = quantile(ili_hat - 0.049, 0.25, na.rm = TRUE),
    q0.75 = quantile(ili_hat - 0.049, 0.75, na.rm = TRUE)
  )

sample_summaries <- sample_summaries %>%
  mutate(
    region = location
  )

us_flu <- us_flu  %>%
  mutate(
    region = ifelse(
      region == "National",
      "US National",
      paste0("HHS ", region)
    )
  )

ggplot() +
  geom_line(
    data = us_flu %>%
      filter(season != predict_season, !is.na(weighted_ili)),
    mapping = aes(x = season_week, y = weighted_ili, group = season), color = "cornflowerblue", alpha = 0.4) +
  geom_line(data = sample_summaries, mapping = aes(x = season_week, y = median_ili_hat)) +
  geom_point(data = sample_summaries, mapping = aes(x = season_week, y = median_ili_hat)) +
  geom_ribbon(data = sample_summaries, mapping = aes(x = season_week, ymin = q0.005, ymax = q0.995), alpha = 0.1) +
  geom_ribbon(data = sample_summaries, mapping = aes(x = season_week, ymin = q0.025, ymax = q0.975), alpha = 0.1) +
  geom_ribbon(data = sample_summaries, mapping = aes(x = season_week, ymin = q0.10, ymax = q0.90), alpha = 0.1) +
  geom_ribbon(data = sample_summaries, mapping = aes(x = season_week, ymin = q0.25, ymax = q0.75), alpha = 0.1) +
  geom_line(
    data = us_flu %>% filter(season == predict_season, !is.na(weighted_ili)),
    mapping = aes(x = season_week, y = weighted_ili), color = "orange", linetype = 2) +
  facet_wrap( ~ region) +
  theme_bw()

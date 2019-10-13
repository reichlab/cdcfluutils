#' Regional influenza incidence in the US (1997 - 2018)
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with 12,056 observations on weighted influenza-like illness 
#' measurements from all HHS regions, including the national level.
#' @source The cdcfluview R package. 
#' @docType data
#' @name flu_data
#' @usage data(flu_data)
NULL

#' Regional influenza incidence in the US (1997 - 2019) including all reports when available.
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with observations on weighted influenza-like illness 
#' measurements from all HHS regions, including the national level.
#' @source The cdcfluview R package. 
#' @docType data
#' @name flu_data_with_backfill
#' @usage data(flu_data_with_backfill)
NULL

#' Flu season "onset thresholds" from the US CDC.
#' 
#' @format A data.frame with 220 observations of the seasonal baseline threshold used for 
#' determining the "flu season onset" in each region of the US.
#' @source \url{https://github.com/cdcepi/FluSight-forecasts/blob/master/wILI_Baseline.csv}
#' @docType data
#' @name flu_onset_baselines
#' @usage data(flu_onset_baselines)
NULL

#' State-level influenza incidence in the US (1997 - 2018)
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with 22,310 observations on influenza-like illness 
#' measurements from all US states.
#' @source The cdcfluview R package. 
#' @docType data
#' @name state_flu_data
#' @usage data(state_flu_data)
NULL
#' Regional influenza incidence in the US (1997 - 2019)
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with 12,650 observations on weighted influenza-like illness 
#' measurements from all HHS regions, including the national level.
#' @source The cdcfluview R package. 
#' @docType data
#' @name flu_data
#' @usage data(flu_data)
NULL

#' National and Regional influenza incidence in the US (1997 - 2019) including all reports when available.
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with observations on weighted influenza-like illness 
#' measurements from all HHS regions and the national level.
#' @source The cdcfluview R package. 
#' @docType data
#' @name nat_reg_flu_data_with_backfill
#' @usage data(nat_reg_flu_data_with_backfill)
NULL

#' State and Local influenza incidence in the US (1997 - 2019) including all reports when available.
#' 
#' A dataset of public influenza data from the US CDC.
#' 
#' @format A data.frame with observations on weighted influenza-like illness 
#' measurements from all states, some territories, and some cities
#' @source The cdcfluview R package. 
#' @docType data
#' @name nat_reg_flu_data_with_backfill
#' @usage data(nat_reg_flu_data_with_backfill)
NULL

#' Flu season "onset thresholds" from the US CDC.
#' 
#' @format A data.frame with 220 observations of the seasonal baseline threshold used for 
#' determining the "flu season onset" in each region of the US.  Note that thresholds
#' for the 1997/1998 through 2006/2007 seasons have been imputed as the mean of later seasons.
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

#' A list of historical revisions to ilinet data at the national and regional
#' level in the US (2003 - 2019).  Entries are revisions to incidence observed
#' so far in the season, treating epidemic week 1 as the first week of the
#' season.
#' 
#' @format list of length 36, where each component is a named list with:
#'   - revision_length: integer. how long the sequence of revisions is
#'   - region_years: data frame. for each revision, what region and year did it
#'     come from
#'   - deltas: numeric matrix with number of rows equal to number of rows of
#'     region_years and number of columns equal to revision_length.
#' @source Postprocessing from the Epidata API
#' @docType data
#' @name ILI_revision_tables
#' @usage data(ILI_revision_tables)
NULL

#' Variances of historical revisions to ilinet data
#' 
#' @format matrix
#' @source Postprocessing from the Epidata API
#' @docType data
#' @name historical_vars_regional
#' @usage data(historical_vars_regional)
NULL
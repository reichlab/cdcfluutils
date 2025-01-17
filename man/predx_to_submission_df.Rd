% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predx-utils.R
\name{predx_to_submission_df}
\alias{predx_to_submission_df}
\title{Converts predx format forecasts into a data frame suitable for writing out
a csv file for submission}
\usage{
predx_to_submission_df(predx_df, ew, year, team = "Kernel of Truth")
}
\arguments{
\item{predx_df}{predx data frame object with samples for submission}

\item{ew}{epidemic week, character in format "01"}

\item{year}{year, numeric or character}

\item{team}{team name, character}
}
\value{
path to csv file
}
\description{
Converts predx format forecasts into a data frame suitable for writing out
a csv file for submission
}

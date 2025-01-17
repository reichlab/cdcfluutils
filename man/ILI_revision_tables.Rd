% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ILI_revision_tables}
\alias{ILI_revision_tables}
\title{A list of historical revisions to ilinet data at the national and regional
level in the US (2003 - 2019).  Entries are revisions to incidence observed
so far in the season, treating epidemic week 1 as the first week of the
season.}
\format{list of length 36, where each component is a named list with:
  - revision_length: integer. how long the sequence of revisions is
  - region_years: data frame. for each revision, what region and year did it
    come from
  - deltas: numeric matrix with number of rows equal to number of rows of
    region_years and number of columns equal to revision_length.}
\source{
Postprocessing from the Epidata API
}
\usage{
data(ILI_revision_tables)
}
\description{
A list of historical revisions to ilinet data at the national and regional
level in the US (2003 - 2019).  Entries are revisions to incidence observed
so far in the season, treating epidemic week 1 as the first week of the
season.
}

library(dplyr)
library(hforecast)

#us_nat_flu <- hforecast::assemble_nat_data() %>%
us_nat_flu <- hforecast::assemble_state_hhs_data(
  region = "national", subset_to_target_weeks = FALSE, include_strains = TRUE
)

us_reg_flu <- hforecast::assemble_state_hhs_data(
  region = "hhs", subset_to_target_weeks = FALSE, include_strains = TRUE
)

cols_both <- colnames(us_nat_flu)[colnames(us_nat_flu) %in% colnames(us_reg_flu)]

us_flu <- us_nat_flu %>%
  dplyr::select(cols_both) %>%
  dplyr::bind_rows(
    us_reg_flu %>%
      dplyr::select(cols_both)
  )

saveRDS(us_flu, file = "inst/make-submission-files/make-GPthingy-submission-file/us_flu.rds")

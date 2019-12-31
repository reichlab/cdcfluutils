library(cdcfluutils)
library(gridExtra)
data <- download_and_preprocess_flu_data()

for(backfill_method in c("none", "post_hoc")) {
  method <- paste0("kcde_backfill_", backfill_method)
  submissions_save_path <- paste0("inst/submissions/region-", method)
  res_file <- paste0(
    submissions_save_path,
    "/csv/",
    "EW", sprintf("%02d", tail(data$week, 1)),
    "-", tail(data$year, 1),
    "-ReichLab_", method,
    ".csv")
  plot_save_path <- paste0(
    submissions_save_path,
    "/plots/",
    "EW", sprintf("%02d", tail(data$week, 1)),
    "-", tail(data$year, 1),
    "-ReichLab_", method,
    ".pdf")
  
  make_predictions_plots(
    preds_save_file = res_file,
    plots_save_file = plot_save_path,
    data = data
  )
}


for(seasonal_difference in c("FALSE", "TRUE")) {
  for(backfill_method in c("none", "post_hoc")) {
    method <- paste0("sarima_seasonal_difference_", seasonal_difference,
           "_backfill_", gsub("-", "_", backfill_method))
    submissions_save_path <- paste0("inst/submissions/region-", method)
    
    res_file <- paste0(
      submissions_save_path,
      "/csv/",
      "EW", sprintf("%02d", tail(data$week, 1)),
      "-", tail(data$year, 1),
      "-ReichLab_", method,
      ".csv")
    plot_save_path <- paste0(
      submissions_save_path,
      "/plots/",
      "EW", sprintf("%02d", tail(data$week, 1)),
      "-", tail(data$year, 1),
      "-ReichLab_", method,
      ".pdf")
    
    make_predictions_plots(
      preds_save_file = res_file,
      plots_save_file = plot_save_path,
      data = data
    )
  }
}

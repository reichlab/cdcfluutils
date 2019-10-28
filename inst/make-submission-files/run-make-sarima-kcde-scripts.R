library(doParallel)

cmds <- c()
for(backfill_method in c("none", "post-hoc")) {
  for(seasonal_difference in c("FALSE", "TRUE")) {
    cmd <- paste0(
      "R CMD BATCH --vanilla \'--args ",
      seasonal_difference, " ",
      backfill_method, " 1",
      "\' inst/make-submission-files/make-regional-sarima-submission-files.R ",
      "inst/make-submission-files/output/output-sarima-test-prediction-step-", seasonal_difference, backfill_method, ".Rout")
    cmds <- c(cmds, cmd)
  }
  
  cmd <- paste0(
    "R CMD BATCH --vanilla \'--args ",
    backfill_method,
    "\' inst/make-submission-files/make-regional-kcde-submission-files.R ",
    "inst/make-submission-files/output/output-kcde-test-prediction-step-", backfill_method, ".Rout")
  cmds <- c(cmds, cmd)
}

registerDoParallel(6)

foreach(i = seq_along(cmds)) %dopar% {
  system(cmds[i])
}

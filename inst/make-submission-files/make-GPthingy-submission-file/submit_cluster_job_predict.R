# R code for SIR fit
library(dplyr)
library(hforecast)

us_flu <- readRDS("inst/make-submission-files/make-GPthingy-submission-file/us_flu.rds")

season <- "2019/2020"
season_week <- us_flu %>% filter(season == UQ(season)) %>% pull(season_week) %>% max()

cores_req <- "1"
mem_req <- "5000"
time_req <- "200:00"
queue_req <- "long"

num_chains <- 50
tree_depth <- 12

save_path <- "/home/er71a/cdcfluutils/inst/make-submission-files/make-GPthingy-submission-file/"
output_path <- "/home/er71a/cdcfluutils/inst/make-submission-files/make-GPthingy-submission-file/cluster-output"
lsfoutfilename <- "predict-GPthingy.out"

#for(chain_num in seq_len(num_chains)) {
for(chain_num in seq_len(2)) {
  case_str <- paste0(gsub("/", "_", season), "_", season_week, "_chain", chain_num)
  filename <- paste0(output_path, "/submit-GPthingy", case_str, ".sh")

  requestCmds <- "#!/bin/bash\n"
  requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
	  "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
	  "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
	  "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
	  "#BSUB -W ", time_req, " # run time\n",
	  "#BSUB -q ", queue_req, " # which queue we want to run in\n")

  cat(requestCmds, file = filename)
  cat("module load gcc/8.1.0\n", file = filename, append = TRUE)
  cat("module load R/3.5.2_gcc8.1.0\n", file = filename, append = TRUE)
  cat(paste0("R CMD BATCH --vanilla \'--args ",
      season, " ",
      season_week, " ",
      tree_depth, " ",
      chain_num,
      "\' /home/er71a/cdcfluutils/inst/make-submission-files/make-GPthingy-submission-file/predict_GPthingy.R ",
      output_path, "/output-GPthingy-", case_str, ".Rout"),
      file = filename, append = TRUE)
  
  bsubCmd <- paste0("bsub < ", filename)
  system(bsubCmd)
}

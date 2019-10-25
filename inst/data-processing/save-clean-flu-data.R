## download and clean regional US flu data
## October 5 2016 - Nicholas Reich - created
## October 7 2016 - Evan Ray - calculate time using MMWRweek package, fix bug in
##   computation of season week
## October 7 2016 - Nicholas Reich - merged US and Regional data into one file.
## October 9 2018 - Nicholas Reich - migrated over to the functionalized form

library(cdcfluview)

flu_data <- download_and_preprocess_flu_data(latest_year=2019)

write.csv(flu_data, file = "data-raw/flu_data.csv")
save(flu_data, file = "data/flu_data.rdata")

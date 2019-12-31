## create long flu baseline dataset
## Nicholas Reich
## October 2016 (updated Oct 2019)

library(tidyverse)

dat <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecasts/master/wILI_Baseline.csv")
dat <- as_tibble(dat) %>%
    gather(key=year, value=baseline, -X) %>%
    transmute(region = factor(as.character(X), levels=c("National", paste0("Region", 1:10))),
              season = str_replace(substr(year, start=2, stop=10), "\\.", "/"), ## removes X at beginning of season name and changes . to /
              baseline=baseline)

## fill in baselines for seasons before 2007/2008 with means
mean_region_baseline <- dat %>%
  group_by(region) %>%
  summarize(baseline = round(mean(baseline), 1))

for(first_season_year in 2006:1997) {
  dat <- bind_rows(
    data.frame(
      season = paste0(first_season_year, "/", first_season_year + 1),
      mean_region_baseline
    ),
    dat
  )
}

flu_onset_baselines <- as.data.frame(dat)

# ggplot(flu_onset_baselines) +
#   geom_line(aes(x=as.numeric(substr(season, 1,4)), y=baseline, color=region))

write.csv(dat, "data-raw/flu_onset_baselines_with_imputations.csv", quote = FALSE, row.names = FALSE)
save(flu_onset_baselines, file = "data/flu_onset_baselines.rdata")

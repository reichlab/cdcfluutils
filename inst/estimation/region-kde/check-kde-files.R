## code to check 2018/2019 submissions

library(FluSight)
library(ggplot2)
library(gridExtra)

all_kde_files <- list.files("inst/submissions/region-kde", full.names = TRUE)

pdf("inst/estimation/region-kde/check_all_files.pdf", width = 12)
for(filename in all_kde_files) {
  print(filename)
  tmp <- read_entry(filename)
  verify_entry(tmp)
  for(reg in unique(tmp$location)){
    p_onset <- plot_onset(tmp, region = reg) + ylim(0,1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size=5))
    p_peakpct <- plot_peakper(tmp, region = reg) + ylim(0,1)
    p_peakwk <- plot_peakweek(tmp, region = reg) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size=5))
    p_1wk <- plot_weekahead(tmp, region = reg, wk = 1, ilimax=13, years = 2017, plot_current=FALSE) + 
      ggtitle(paste(reg, ": 1 wk ahead")) + ylim(0,1)
    p_2wk <- plot_weekahead(tmp, region = reg, wk = 2, ilimax=13, years = 2017, plot_current=FALSE) + ylim(0,1)
    p_3wk <- plot_weekahead(tmp, region = reg, wk = 3, ilimax=13, years = 2017, plot_current=FALSE) + ylim(0,1)
    p_4wk <- plot_weekahead(tmp, region = reg, wk = 4, ilimax=13, years = 2017, plot_current=FALSE) + ylim(0,1)
    grid.arrange(p_1wk, p_2wk, p_3wk, p_4wk, p_onset, p_peakpct, p_peakwk, ncol=4)
  }
}
dev.off()
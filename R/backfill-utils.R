investigate_backfill <- function(){
  library(dplyr)
  library(ggplot2)
  library(forecast)
  library(cowplot)
  ### READ DATA AND PLOT
  backfill_data <- readRDS("data/flu_data_with_backfill.rds")
  # SUBSET NATIONAL DATA
  backfill_data <- backfill_data[backfill_data$region=="nat",]
  backfill_data$release_date <- as.Date(backfill_data$release_date)
  backfill_data <- backfill_data[order(backfill_data$epiweek),]
  
  #PLOT FULL DATA
  
  ## CREATE FULLY OBSERVED
  fully_observed_data <-backfill_data %>% group_by(region,epiweek) %>% summarize(wili = last(wili) )
  plot(fully_observed_data$wili,type='l')
  ## CREATE PARTIALLY OBSERVED
  
  partially_observed_data <-backfill_data %>% group_by(release_date,region,epiweek) %>% summarize(wili = last(wili) )
  #DEFINE A SEQUENCE OF TEST YEARS
  years<- seq(2003,2005,by=1)
  delta_trajectories <- matrix(NA,nrow=length(years)*35,ncol=35)
  row_iter <-1
  for (j in 1:(length(years)-1)){
    plot_epi_week <- as.numeric(paste(years[j+1],"19",sep=""))
    partially_observed_week_dates <- seq(as.Date(paste(years[j],"-12-01",sep="")),as.Date(paste(years[j+1],"-6-01",sep="")),by=7)
    fully_observed_wili <- fully_observed_data[fully_observed_data$epiweek<=plot_epi_week & fully_observed_data$epiweek >= as.numeric(paste(years[j],"40",sep="")),]$wili
    for(i  in 1:length(partially_observed_week_dates)){
      partially_observed_week <- partially_observed_week_dates[i]
      tmp <- partially_observed_data[partially_observed_data$epiweek >= as.numeric(paste(years[j],"40",sep="")) & partially_observed_data$epiweek<=plot_epi_week & partially_observed_data$release_date <= partially_observed_week,]
      tmp <- tmp  %>% group_by(region,epiweek) %>% summarize(wili = wili[which.max(release_date)])
      p<- ggplot(tmp,aes(x=1:nrow(tmp),y=wili-fully_observed_wili[1:length(tmp$wili)]),col='red') + ylab("Difference from truth and currently reported wILI") + xlab(paste("Time from",years[j],"40")) + geom_line() + ylim(-1.5,1.5) + xlim(1,35)+ ggtitle(years[j]) +theme_bw()
      #p <- p + geom_vline(xintercept=7,alpha=.4,linetype="dotted",col="red")
      #ggsave(paste(years[j],"-",i,".pdf",sep=""),p,scale = 1.5)
      tmp_trajectory <- tmp$wili-fully_observed_wili[1:length(tmp$wili)]
      while (length(tmp_trajectory) <35){
        tmp_trajectory <- c(tmp_trajectory,NA)
      }
      delta_trajectories[row_iter,] <-tmp_trajectory
      row_iter <- row_iter +1
    }
  }
  
}



#' Create non-parametric matrix of backfill with dimensions
#'  n_seasons x n_weeks x n_lags
#'

create_e_matrix <- function(){
  
  ggplot(weekILI[weekILI$epiweek < 200352,],aes(x=lag,y=wili,col=factor(epiweek))) + geom_line() + xlab("Lag in Season Week")
  backfill_data <- readRDS("data/flu_data_with_backfill.rds")
  
}

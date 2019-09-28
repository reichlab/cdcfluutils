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



library(plyr) # for rbind.fill
library(dplyr)
library(tidyr)
library(MMWRweek)
library(jsonlite)
source("https://raw.githubusercontent.com/cmu-delphi/delphi-epidata/master/src/client/delphi_epidata.R")
move_k_week_ahead <- function(epiweek,k){
  ### UTILITY FUNCTION TO MOVE K WEEK AHEAD
  if (k==1){
    if (substr(epiweek,5,7) == 52){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (as.numeric(substr(epiweek,5,7)) < 9){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+1)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 1
      new_year <- substr(epiweek,1,4)
    }
  } else if (k==2){
    if (substr(epiweek,5,7) == 51){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 52){
      new_week <- "02"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    }else if (as.numeric(substr(epiweek,5,7)) < 8){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+2)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 2
      new_year <- substr(epiweek,1,4)
    }
  } else if (k==3){
    if (substr(epiweek,5,7) == 50){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 51){
      new_week <- "02"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 52){
      new_week <- "03"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (as.numeric(substr(epiweek,5,7)) < 7){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+3)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 3
      new_year <- substr(epiweek,1,4)
    }
  } else if (k==4){
    if (substr(epiweek,5,7) == 49){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 50){
      new_week <- "02"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 51){
      new_week <- "03"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 52){
      new_week <- "04"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (as.numeric(substr(epiweek,5,7)) < 6){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+4)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 4
      new_year <- substr(epiweek,1,4)
    }
  }
  return (paste0(new_year,new_week))
}


rRevisedILI <- function(n, observed_inc, epiweek_idx, region, season) {
  #' sample revised ILI trajectories
  #'
  #' @param n number of revised trajectories to sample
  #' @param observed_inc observed incidence so far this season, starting at EW40 and going to the most recent report
  #' @param epiweek most recent epidemic week (equals 40 + length(observed_inc) - # weeks in season)
  #' @param region region trajectory was made from 
  #' @param season current season in 20xx/20xx+1 format
  #' @return n by length(observed_inc) matrix of samples for possible revised ili.
  
  ##delayed_data 
  flu_data_with_backfill <- readRDS("data/flu_data_with_backfill.rds")
  flu_data_with_backfill$issue_week <- unlist(lapply(flu_data_with_backfill$issue,function(x){return(substr(x,5,7))}))
  flu_data_with_backfill$season <- unlist(lapply(flu_data_with_backfill$issue,function(x){
    
    current_epiweek <- as.numeric(substr(x,5,7))
    if (current_epiweek <= 20){
      return(as.numeric(substr(x,1,4))-1)
    } else{
      return(as.numeric(substr(x,1,4)))
    }
    
  }))
  
  if (epiweek_idx <= 20){
    time_in <- 12 + epiweek_idx
  } else{
    time_in <- epiweek_idx -40 + 1
  }
  # fully observed data
  regions <- c(paste0("hhs",1:10),"nat")
  fully_observed_data <- flu_data_with_backfill %>% group_by(epiweek,region) %>% filter(lag==max(lag))  %>% arrange(epiweek)
  total_traj <-  matrix(nrow = n,ncol= time_in)
  for (i in 1:n){
    sampled_region <- sample(regions,1)
    sampled_season <- sample(paste0(20,10:as.numeric(substr(season,2,4))),1)
    sampled_epiweek <- paste0(sampled_season,epiweek_idx)
    avail <- flu_data_with_backfill[flu_data_with_backfill$region == sampled_region & flu_data_with_backfill$issue <= sampled_epiweek,] %>% group_by(epiweek) %>% filter(lag ==max(lag)) %>% arrange(epiweek)
    fully_obs <- fully_observed_data[fully_observed_data$region == sampled_region & fully_observed_data$epiweek <= sampled_epiweek,]
    delta <- tail(avail$wili,time_in)-tail(fully_obs$wili,time_in)
    total_traj[i,] <- observed_inc-delta
  }
  
  ## add nowcast 
  
  current_season_epiweek <- ifelse(epiweek_idx <=20,paste0(substr(season,6,12),epiweek_idx),paste0(substr(season,1,4),epiweek_idx))
  
  req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,1)))
  nowcast_json <- jsonlite::prettify(rawToChar(req$content))
  nowcast_obj_1_wk_ahead <- fromJSON(nowcast_json)
  
  req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,2)))
  nowcast_json <- jsonlite::prettify(rawToChar(req$content))
  nowcast_obj_2_wk_ahead <- fromJSON(nowcast_json)
  
  total_traj <- cbind(total_traj, rep(nowcast_obj_1_wk_ahead$epidata$value,n),rep(nowcast_obj_2_wk_ahead$epidata$value,n))
  return (total_traj)
  
  
}


#' Utility function to move k week ahead
move_k_week_ahead <- function(epiweek,k){
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


#' sample revised ILI trajectories
#'
#' @param n number of revised trajectories to sample
#' @param observed_inc observed incidence so far this season, starting at EW40 and going to the most recent report
#' @param epiweek most recent epidemic week (equals 40 + length(observed_inc) - # weeks in season)
#' @param region region trajectory was made from 
#' @param season current season in 20xx/20xx+1 format
#' @param add_nowcast logical; add nowcast based on delphi epicast?
#' @return n by length(observed_inc) matrix of samples for possible revised ili.
#' 
#' @export
rRevisedILI <- function(n, observed_inc, epiweek_idx, region, season, add_nowcast = FALSE) {
  ##delayed_data 
  flu_data_with_backfill$issue_week <- unlist(lapply(
    flu_data_with_backfill$issue,
    function(x) { return(substr(x,5,7)) }
  ))
  flu_data_with_backfill$season <- unlist(lapply(
    flu_data_with_backfill$issue,
    function(x) {
      current_epiweek <- as.numeric(substr(x,5,7))
      if (current_epiweek <= 20){
        return(as.numeric(substr(x,1,4))-1)
      } else{
        return(as.numeric(substr(x,1,4)))
      }
    }
  ))
  
  if (epiweek_idx <= 20){
    time_in <- 12 + epiweek_idx
  } else{
    time_in <- epiweek_idx - 40 + 1
  }
  
  # fully observed data
  regions <- c(paste0("hhs",1:10),"nat")
  fully_observed_data <- flu_data_with_backfill %>%
    dplyr::group_by(epiweek, region) %>%
    dplyr::filter(lag == max(lag)) %>%
    dplyr::arrange(epiweek)
  
  total_traj <-  matrix(nrow = n, ncol= time_in)
  for (i in 1:n) {
    sampled_region <- sample(regions,1)
    sampled_season <- sample(paste0(20,10:as.numeric(substr(season,2,4))),1)
    sampled_epiweek <- paste0(sampled_season,epiweek_idx)
    avail <- flu_data_with_backfill[flu_data_with_backfill$region == sampled_region & flu_data_with_backfill$issue <= sampled_epiweek,] %>% group_by(epiweek) %>% filter(lag ==max(lag)) %>% arrange(epiweek)
    fully_obs <- fully_observed_data[fully_observed_data$region == sampled_region & fully_observed_data$epiweek <= sampled_epiweek,]
    delta <- tail(avail$wili,time_in)-tail(fully_obs$wili,time_in)
    total_traj[i,] <- observed_inc-delta
  }
  
  ## add nowcast 
  if(add_nowcast) {
    current_season_epiweek <- ifelse(epiweek_idx <=20,paste0(substr(season,6,12),epiweek_idx),paste0(substr(season,1,4),epiweek_idx))
    
    req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,1)))
    nowcast_json <- jsonlite::prettify(rawToChar(req$content))
    nowcast_obj_1_wk_ahead <- fromJSON(nowcast_json)
    
    req <- curl::curl_fetch_memory(url=paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=nowcast&locations=",region,"&epiweeks=",move_k_week_ahead(current_season_epiweek,2)))
    nowcast_json <- jsonlite::prettify(rawToChar(req$content))
    nowcast_obj_2_wk_ahead <- fromJSON(nowcast_json)
    
    total_traj <- cbind(total_traj, rep(nowcast_obj_1_wk_ahead$epidata$value,n),rep(nowcast_obj_2_wk_ahead$epidata$value,n))
  }
  
  return (total_traj)
}

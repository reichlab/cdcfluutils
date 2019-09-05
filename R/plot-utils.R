#' Plots Week Ahead Forecasts
#'
#' This function allows you to plot Weak Ahead Predictions
#' @param dat Expects a data csv in the form of a CDC fluview submission see \code{FluSight} package for a minimal submission
#' @param region Specifies the region to be plotted
#' @param wk Numeric. How many weeks ahead to plot. Defaults to 1.
#' @param ilimax Numeric. Max level of ILI percentage to be plotted
#' @keywords Week Ahead Prediction Plots
#' @export
#' @examples
#' plotWeekAhead()
my_plot_weekahead <- function(dat, region, wk=1, ilimax, years=NA){
    require(ggplot2)
    require(cdcfluview)
    
    d <- suppressWarnings(subset(dat, location==region & as.numeric(as.character(bin_start_incl)) <= ilimax))
    d <- d[grep("wk ahead", d$target),]
    
    NatFluDat <- get_flu_data("national", NA, "ilinet", years=years)
    RegFluDat <- get_flu_data("hhs", 1:10, "ilinet", years=years)
    StateFluDat <- get_flu_data("state", "all", "ilinet", years=years)
    ## kludgy work-around
    StateFluDat$`% WEIGHTED ILI` <- StateFluDat$`%UNWEIGHTED ILI`
    NatFluDat$REGION[NatFluDat$REGION=="X"] <- "National"
    CFluDat <- rbind(
        tail(RegFluDat, n=10),
        tail(NatFluDat, n=1),
        tail(StateFluDat, n=53))[,c("REGION","% WEIGHTED ILI")]
    if(region=="US National"){
        CurrentILIPer <- as.numeric(CFluDat[CFluDat$REGION=="National", "% WEIGHTED ILI"])
    }else if(region %in% paste("Region", 1:10)){
        CurrentILIPer <- as.numeric(CFluDat[CFluDat$REGION==paste("Region", strsplit(region, " ")[[1]][3]), "% WEIGHTED ILI"])
    } else {
        CurrentILIPer <- as.numeric(CFluDat[CFluDat$REGION==region, "% WEIGHTED ILI"])
    }
    
    if(wk==1){
        return(ggplot(data=subset(d, target=="1 wk ahead"), aes(x=as.numeric(as.character(bin_start_incl)), y=value)) +
                geom_point() + labs(title = "1 Week Ahead", x="ILI%", y="Prob") + geom_vline(aes(xintercept = CurrentILIPer)))
    }
    
    
    if(wk==2){
        return(ggplot(data=subset(d, target=="2 wk ahead"), aes(x=as.numeric(as.character(bin_start_incl)), y=value)) +
                geom_point() + labs(title = "2 Week Ahead", x="ILI%", y="Prob") + geom_vline(aes(xintercept = CurrentILIPer)))
    }
    
    if(wk==3){
        return(ggplot(data=subset(d, target=="3 wk ahead"), aes(x=as.numeric(as.character(bin_start_incl)), y=value)) +
                geom_point() + labs(title = "3 Week Ahead", x="ILI%", y="Prob") + geom_vline(aes(xintercept = CurrentILIPer)))
    }
    
    if(wk==4){
        return(ggplot(data=subset(d, target=="4 wk ahead"), aes(x=as.numeric(as.character(bin_start_incl)), y=value)) +
                geom_point() + labs(title = "4 Week Ahead", x="ILI%", y="Prob") + geom_vline(aes(xintercept = CurrentILIPer)))
    }
    
}

#' Make plots of prediction submissions for flu contest: so far, seasonal targets only
#'
#' @param preds_save_file path to a file with predictions in csv submission format
#' @param plots_save_file path to a pdf file where plots should go
#' @param data data observed so far this season
#' @export
make_predictions_plots <- function(
  preds_save_file,
  plots_save_file,
  data
) {
  require("grid")
  require("ggplot2")
  
  predictions <- read.csv(preds_save_file)
  regional <- data$region_type[1] =="HHS Regions"
  
  if(regional) {
    preds_region_map <- data.frame(
      internal_region = c("National", paste0("Region ", 1:10)),
      preds_region = c("US National", paste0("HHS Region ", 1:10))
    )
  } else {
    preds_region_map <- data.frame(
      internal_region = unique(predictions$Location),
      preds_region = unique(predictions$Location)
    )
  }
  
  current_season <- tail(data$season, 1)
  current_year <- tail(data$year, 1)
  
  pdf(plots_save_file)
  
  for(region in unique(data$region)) {
    preds_region <- preds_region_map$preds_region[preds_region_map$internal_region == region]
    
    ## Observed incidence
    p_obs <- ggplot(data[data$region == region & data$season == current_season, ]) +
      expand_limits(x = c(0, 42), y = c(0, 13)) +
      geom_line(aes(x = season_week, y = weighted_ili)) +
      ggtitle("Observed incidence") +
      theme_bw()
    
    if(regional){
      p_obs <- p_obs +
        geom_hline(yintercept = get_onset_baseline(region, season = current_season), colour = "red")
    }
    
    ## Onset
    if(regional) {
      reduced_preds <- predictions[predictions$Location == preds_region & predictions$Target == "Season onset" & predictions$Type == "Bin", ] %>%
        mutate(
          season_week = year_week_to_season_week(as.numeric(as.character(Bin_start_incl)), current_year)
        )
      point_pred <- predictions[predictions$Location == preds_region & predictions$Target == "Season onset" & predictions$Type == "Point", , drop = FALSE] %>%
        mutate(
          season_week = year_week_to_season_week(Value, current_year)
        )
      p_onset <- ggplot(reduced_preds) +
        geom_line(aes(x = season_week, y = Value)) +
        geom_vline(xintercept = point_pred$season_week, colour = "red") +
        expand_limits(x = c(0, 42)) +
        ylab("predicted probability of onset") +
        ggtitle("Onset") +
        theme_bw()
    }
    
    ## Peak Timing
    reduced_preds <- predictions[predictions$Location == preds_region & predictions$Target == "Season peak week" & predictions$Type == "Bin", ] %>%
      mutate(
        season_week = year_week_to_season_week(as.numeric(as.character(Bin_start_incl)), current_year)
      )
    point_pred <- predictions[predictions$Location == preds_region & predictions$Target == "Season peak week" & predictions$Type == "Point", , drop = FALSE] %>%
      mutate(
        season_week = year_week_to_season_week(Value, current_year)
      )
    p_peak_timing <- ggplot(reduced_preds) +
      geom_line(aes(x = season_week, y = Value)) +
      geom_vline(xintercept = point_pred$season_week, colour = "red") +
      expand_limits(x = c(0, 42)) +
      ylab("predicted probability of peak") +
      ggtitle("Peak timing") +
      theme_bw()
    
    
    ## Peak Incidence
    reduced_preds <- predictions[predictions$Location == preds_region & predictions$Target == "Season peak percentage" & predictions$Type == "Bin", ] %>%
      mutate(inc_bin = as.numeric(as.character(Bin_start_incl)))
    point_pred <- predictions[predictions$Location == preds_region & predictions$Target == "Season peak percentage" & predictions$Type == "Point", , drop = FALSE] %>%
      mutate(inc_bin = Value)
    p_peak_inc <- ggplot(reduced_preds) +
      geom_line(aes(x = inc_bin, y = Value)) +
      geom_vline(xintercept = point_pred$inc_bin, colour = "red", data) +
      expand_limits(x = c(0, 13)) +
      ylab("predicted probability of peak incidence") +
      coord_flip() +
      ggtitle("Peak incidence") +
      theme_bw()
    
    grid.newpage()
    pushViewport(viewport(layout =
                            grid.layout(nrow = 4,# ifelse(regional, 4, 3), ## adjustment for onset
                                        ncol = 2,
                                        heights = unit(c(2, 1, 1, 1), c("lines", "null", "null", "null")))))
    
    grid.text(preds_region,
              gp = gpar(fontsize = 20),
              vp = viewport(layout.pos.col = 1:2, layout.pos.row = 1))
    if(regional){
      print(p_onset, vp = viewport(layout.pos.col = 1, layout.pos.row = 2))
    }
    print(p_obs, vp = viewport(layout.pos.col = 1, layout.pos.row = 3))
    print(p_peak_timing, vp = viewport(layout.pos.col = 1, layout.pos.row = 4))
    print(p_peak_inc, vp = viewport(layout.pos.col = 2, layout.pos.row = 3))
    
    recent_obs <- data[data$region == region & data$season == current_season, "weighted_ili"]
    recent_obs <- tail(recent_obs, 1)
    p_1wk <- my_plot_weekahead(res, region = preds_region, wk = 1, ilimax=13, years = 2018:2019) + ggtitle(paste(preds_region, ": 1 wk ahead")) + ylim(0,1) + geom_vline(xintercept = recent_obs)
    p_2wk <- my_plot_weekahead(res, region = preds_region, wk = 2, ilimax=13, years = 2018:2019) + ylim(0,1) + geom_vline(xintercept = recent_obs)
    p_3wk <- my_plot_weekahead(res, region = preds_region, wk = 3, ilimax=13, years = 2018:2019) + ylim(0,1) + geom_vline(xintercept = recent_obs)
    p_4wk <- my_plot_weekahead(res, region = preds_region, wk = 4, ilimax=13, years = 2018:2019) + ylim(0,1) + geom_vline(xintercept = recent_obs)
    grid.arrange(p_1wk, p_2wk, p_3wk, p_4wk, ncol=1)
  }
  
  dev.off()
}

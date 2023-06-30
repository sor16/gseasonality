#' Summary method for seasonality models
#'
#' Summarize a seasonality gam model object
#' @param object an object of class "seasm".
#' @param ... Not used for this function
#' @seealso \code{\link{seasonality_gam}},
#' @export
summary.seasm <- function(object,...){
  ptr_summary <- object$summary$ptr
  names(ptr_summary) <- paste0(names(ptr_summary),c('','-2.5%','-97.5%'))
  cat("\nPeak-to-Trough-Ratio (PTR):\n")
  print(format(ptr_summary,digits=2,nsmall=2),row.names=F,right=T)
  cat("\nSeasonal characteristics:\n")
  print(format(object$summary$peak_trough,digits=2,nsmall=2),right=T)
  cat("\nCross-over points:\n")
  cat(paste0('low_to_high ',format(object$summary$cross_over_points[['low_to_high']],digits=2,nsmall=2)))
  cat(paste0('\nhigh_to_low ',format(object$summary$cross_over_points[['high_to_low']],digits=2,nsmall=2)))
  cat("\n\nDeviance explained:",paste0(format(100*object$summary$deviance['deviance_expl'],digits=2,nsmall=2),'%'))
  cat("\nDispersion:",format(object$summary$dispersion,digits=2,nsmall=2))
}

#' Print method for seasonality models
#'
#' Print the summary of the seasonality model
#' @param x an object of class "seasm".
#' @param ... Not used for this function
#' @seealso \code{\link{seasonality_gam}}
#' @export
print.seasm <- function(x,...){
  summary(x)
}

#' @importFrom mgcv gam
#' @importFrom dplyr n bind_rows tibble mutate
#' @importFrom lubridate ymd ym
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_rect geom_hline scale_fill_manual scale_x_date scale_x_continuous scale_y_continuous ggtitle xlab ylab theme theme_bw theme_classic
#' @importFrom rlang .data
plot_fun <- function(object,type,...){
  if(type=='fit'){
    monthly_counts_endpoint <- object$monthly_counts
    years <- unique(monthly_counts_endpoint$EVENT_YEAR)
    season_dat <- bind_rows(list(summer=tibble(start=ym(paste0(years,'-',5)),
                                               end=ym(paste0(years,'-',9))),
                                 winter=tibble(start=ym(paste0(years,'-',10)),
                                               end=ym(paste0(years+1,'-',4)))),.id='season')
    season_dat <- bind_rows(season_dat,tibble(season='winter',
                                              start=ym(paste0(years[1],'-',1)),
                                              end=ym(paste0(years[1],'-',4))))
    monthly_counts_endpoint <- left_join(object$fit_table,object$monthly_counts,by=c('EVENT_YEAR','EVENT_MONTH')) %>%
                                mutate(EVENT_DATE=ymd(paste0(EVENT_YEAR,'-',EVENT_MONTH,'-','01')),
                                       month=1:n())
    ggplot(monthly_counts_endpoint) +
      geom_point(aes(x=.data$EVENT_DATE,
                     y=.data$COUNT),alpha=0.5) +
      geom_line(aes(x=.data$EVENT_DATE,
                    y=exp(.data$fit))) +
      geom_line(aes(x=.data$EVENT_DATE,
                    y=exp(fit_null)),
                col='red') +
      geom_rect(data=season_dat,aes(xmin=start-50,
                                    xmax=end,
                                    ymin=0.97*min(monthly_counts_endpoint$COUNT,na.rm=T),
                                    ymax=1.03*max(monthly_counts_endpoint$COUNT,na.rm=T),
                                    fill=season),alpha=0.1) +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(expand=c(0.005,0),date_breaks='3 years',date_labels='%Y') +
      scale_fill_manual(values=c('gold3','turquoise')) +
      ggtitle('') +
      xlab('') +
      ylab('Count') +
      theme_classic() +
      theme(legend.position='none')
  }else if(type=='seasonality'){
    ggplot(object$seasonality_term) +
      geom_line(aes(x=month,
                    y=est)) +
      geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.3,fill='blue') +
      geom_hline(yintercept=0,linetype='dashed',color='red') +
      scale_x_continuous(breaks=seq(1,12)) +
      xlab('Month number') +
      ylab('Smooth term') +
      ggtitle('Seasonal component') +
      theme_bw()
  }else if(type=='annual'){
    ggplot(object$annual_term) +
      geom_line(aes(x=year,
                    y=est)) +
      geom_ribbon(aes(x=year,ymin=lower,ymax=upper),alpha=0.3,fill='blue') +
      geom_hline(yintercept=0,linetype='dashed',color='red') +
      scale_x_continuous(limits = c(object$year_start,object$year_end),expand=c(0.03,0)) +
      xlab('Year') +
      ylab('Smooth term') +
      ggtitle('Annual component') +
      theme_bw()
  }

}

#' Autoplot method for seasonality GAM on monthly counts
#'
#' Visualize different aspects of the seasonality GAM
#' @param object an object of class "seasm".
#' @param ... other plotting parameters (not used in this function)
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"fit"}{ to plot the GAM fit to the monthly counts}
#'                       \item{"seasonality"}{ to plot the seasonal smooth term in th GAM model}
#'                       \item{"annual"}{ to plot the annual smooth term in th GAM model}
#'                    }
#'
#' @return returns an object of class "ggplot2".
#' @seealso \code{\link{seasonality_gam}} and \code{\link{summary.seasm}}, \code{\link{plot.seasm}}
#' @importFrom ggplot2 autoplot
#' @export
autoplot.seasm <- function(object,...,type='fit'){
  plot_fun(object,type=type)
}

#' Plot method for seasonality GAM on monthly counts
#'
#' Visualize different aspects of the seasonality GAM
#' @param object an object of class "seasm".
#' @param ... other plotting parameters (not used in this function)
#' @param type a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
#'                    \itemize{
#'                       \item{"fit"}{ to plot the GAM fit to the monthly counts}
#'                       \item{"seasonality"}{ to plot the seasonal smooth term in th GAM model}
#'                       \item{"annual"}{ to plot the annual smooth term in th GAM model}
#'                       \item{"panel"}{ to plot all three plots above (default)}
#'                    }
#'
#' @return No return value, called for side effects.
#' @seealso \code{\link{seasonality_gam}}, \code{\link{summary.seasm}}, \code{\link{autoplot.seasm}}
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 autoplot
#' @export
plot.seasm <- function(x,...,type='panel'){
  if(is.null(type) || type!='panel'){
    p <- autoplot(x,type=type)
  }else{
    p <- autoplot(x,type='fit') / (autoplot(x,type='annual') + autoplot(x,type='seasonality'))
  }
  print(p)
}

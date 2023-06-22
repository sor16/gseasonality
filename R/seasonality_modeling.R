#' Run seasonality model using monthly counts data
#'
#' Runs a GAM which includes an annual and a seasonal term assuming a quasipoisson count distribution
#' @param data data.frame with columns:
#'                    \itemize{
#'                       \item{"ID"}{ a character denoting the ID of the patient}
#'                       \item{"EVENT_DATE"}{ The first date when the ICD code is found in the registry}
#'                    }
#' @param year_start integer denoting the starting year of the analysis.
#' Defaults to NULL in which case the starting year is taken to be the earliest year in data
#' @param year_end integer denoting the final year of the analysis.
#' Defaults to NULL in which case the final year is taken to be the latest year in data
#' @param adjusted boolean denoting whether to perform seasonal behavioural adjustment. Defaults to FALSE
#' @param ... Not used for this function
#' @details Details here
#' @return The function returns an object of type "seasm"
#' @seealso \code{\link{summary.seasm} and \link{plot.seasm}} for further studying the output
#' @importFrom dplyr mutate filter count summarise group_by
#' @importFrom lubridate ymd year month
#' @importFrom magrittr "%>%"
#' @export
seasonality_gam <- function(data,year_start=NULL,year_end=NULL,adjusted=F){
  data_filtered <- mutate(data,EVENT_DATE=ymd(EVENT_DATE)) %>%
                   mutate(EVENT_YEAR=year(EVENT_DATE),
                          EVENT_MONTH=month(EVENT_DATE),
                          EVENT_DAY=day(EVENT_DATE))
  if(is.null(year_start)){
    year_start <- min(data_filtered$EVENT_YEAR)
  }
  if(is.null(year_end)){
    year_end <- max(data_filtered$EVENT_YEAR)
  }
  data_filtered <- filter(data_filtered,EVENT_YEAR>=year_start,EVENT_YEAR<=year_end)
  monthly_counts <- count(data_filtered,EVENT_YEAR,EVENT_MONTH,name = 'COUNT')
  mod_list <- run_seasonality_gam(dat=monthly_counts,adjustment=ifelse(adjusted,'binary',''),a=1,b=13)
  seasonality_summary <- summarise_seasonality(mod_list=mod_list,monthly_counts=monthly_counts)
  seasonality_obj <- list()
  attr(seasonality_obj, "class") <- "seasm"
  seasonality_obj$fit <- mod_list$seasonal
  seasonality_obj$fit_null <- mod_list$null
  seasonality_obj$summary <- seasonality_summary
  seasonality_obj$data <- data_filtered
  seasonality_obj$monthly_counts <- monthly_counts
  return(seasonality_obj)
}

#' @importFrom stats approx
get_seasonal_spline_adj <- function(seasonal_spline_avg,months){
  approx(x=seasonal_spline_avg$month,
         y=seasonal_spline_avg$avg_seasonal_val,
         xout=months)$y
}

get_grid <- function(type,year_start,year_end,adjustment='',by_year=0.1,by_month=0.1,seasonal_spline_avg=NULL){
  if(type=='seasonal'){
    grid_dat <-data.frame(EVENT_YEAR=year_start, EVENT_MONTH=seq(0.5,12.5,by=by_month))
  }else if(type=='annual'){
    grid_dat <- data.frame(EVENT_YEAR=seq(year_start,year_end,by=by_year),EVENT_MONTH=1)
  }else{
    stop('Type not recognized')
  }
  grid_dat$nr_days <- 30
  if(grepl('mean',adjustment)){
    grid_dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,EVENT_MONTH)
  }else if(adjustment=='binary'){
    grid_dat$july <- as.integer(floor(EVENT_MONTH)==7)
    grid_dat$december <- as.integer(floor(EVENT_MONTH)==12 | EVENT_MONTH<1)

  }
  return(grid_dat)
}

#' @importFrom mgcv gam
#' @importFrom lubridate ymd year month day
#' @importFrom dplyr mutate tibble count group_by summarise inner_join
#' @importFrom stats quasipoisson as.formula
#' @importFrom magrittr "%>%"
run_seasonality_gam <- function(dat,a,b,adjustment='',adj_month_length=F,mod_null='adj',seasonal_spline_type='cp',k_seasonal=6,seasonal_spline_avg=NULL){
  k_annual=6
  dat$nr_days <- 30
  offset_null <- 'log(nr_days)'
  if(adjustment=='mean_covariate'){
    dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,dat$EVENT_MONTH)
    offset_adj <- ''
    cov_adj <- ' + avg_seasonal_val'
  }else if(adjustment=='binary'){
    dat$july <- as.integer(floor(dat$EVENT_MONTH)==7)
    dat$december <- as.integer(floor(dat$EVENT_MONTH)==12)
    offset_adj <- ''
    cov_adj <- ' + july + december'
  }else if(adjustment=='mean_offset'){
    dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,dat$EVENT_MONTH)
    offset_adj <- ' + avg_seasonal_val'
    cov_adj <- ''
  }else{
    offset_adj <- ''
    cov_adj <- ''
  }
  f_start <- 'COUNT ~ s(EVENT_YEAR, k=k_annual, bs="ps")'
  if(mod_null=='adj'){
    f_null <- paste0(f_start,' + offset(',offset_null,offset_adj,')',cov_adj)
    f_seasonal <- f_null
  }else if(mod_null=='noadj'){
    f_null <- paste0(f_start,' + offset(',offset_null,')')
    f_seasonal <- paste0(f_start,' + offset(',offset_null,offset_adj,')',cov_adj)
  }else{
    stop('mod_null argument not recognized')
  }
  f_seasonal <- paste0(f_seasonal,' + s(EVENT_MONTH, k=k_seasonal, bs="',seasonal_spline_type,'")')
  family <- quasipoisson()
  seasonal <- gam(as.formula(f_seasonal), family=family, data=dat, knots=list(EVENT_MONTH = c(a-0.5, b-0.5)), scale=-1,method = "REML")
  null <- gam(as.formula(f_null),family=family, data=dat, scale=-1,method = "REML")
  return(list(seasonal=seasonal, null=null))
}

#' @importFrom mgcv gam
summarise_seasonality <- function(mod_list,monthly_counts,adjustment='',seasonal_spline_avg=NULL){
  seasonal_grid <- get_grid(type='seasonal',year_start=min(monthly_counts$EVENT_YEAR),year_end=max(monthly_counts$EVENT_YEAR),adjustment=adjustment,seasonal_spline_avg=seasonal_spline_avg,by_month = 0.01)
  seasonal_component_monthly <- predict(mod_list$seasonal,
                                        newdata=seasonal_grid,
                                        type='terms')
  seasonal_component_monthly_est <- seasonal_component_monthly[,'s(EVENT_MONTH)']
  idx_peak <- which.max(seasonal_component_monthly_est)
  val_peak <- exp(seasonal_component_monthly_est[idx_peak])
  idx_trough <- which.min(seasonal_component_monthly_est)
  val_trough <- exp(seasonal_component_monthly_est[idx_trough])
  month_peak <- seasonal_grid$EVENT_MONTH[idx_peak]
  month_trough <- seasonal_grid$EVENT_MONTH[idx_trough]
  month_peak <- ifelse(month_peak<1,month_peak+12,month_peak) # adjust to month ragne [1,13)
  month_trough <- ifelse(month_trough<1,month_trough+12,month_trough) # adjust to month ragne [1,13)
  ptr_est <- val_peak / val_trough
  CI <- get_monte_carlo_PTR_CI(mod=mod_list$seasonal,monthly_counts=monthly_counts,nr_iter=100,adjustment=adjustment,seasonal_spline_avg = seasonal_spline_avg)
  peak_trough <- data.frame(peak=c(month_peak,val_peak),trough=c(month_trough,val_trough))
  rownames(peak_trough) <- c('month','value')
  summary_list <- list(ptr=data.frame(estimate=ptr_est,lower=CI$ptr['2.5%'],upper=CI$ptr['97.5%']),
                       peak_trough=peak_trough,
                       pval=summary(mod_list$seasonal)$s.table['s(EVENT_MONTH)','p-value'],
                       deviance=c('deviance'=mod_list$seasonal$deviance,
                                  'deviance_expl'=(mod_list$null$deviance-mod_list$seasonal$deviance)/mod_list$null$deviance),
                       dispersion=mod_list$seasonal$scale)
  return(summary_list)
}

#' @importFrom mgcv gam
#' @importFrom stats quantile
#' @importFrom MASS mvrnorm
get_monte_carlo_PTR_CI <- function(mod,monthly_counts,nr_iter,adjustment,seasonal_spline_avg){
  relevant <- grep('EVENT_MONTH',names(mod$coefficients))
  V_f <- vcov(mod)[relevant,relevant]
  s_hat <- coef(mod)[relevant]
  prediction_grid <- get_grid(type='seasonal',year_start=min(monthly_counts$EVENT_YEAR),year_end=max(monthly_counts$EVENT_YEAR),by_month = 0.01,adjustment=adjustment,seasonal_spline_avg=seasonal_spline_avg)
  lp <- predict(mod, newdata = prediction_grid , type = "lpmatrix")[,relevant]
  s_sample <- mvrnorm(n=nr_iter,mu=s_hat,Sigma=V_f)
  smooth_sample <- lp %*% t(s_sample)
  month_grid_adjusted <- ifelse(prediction_grid$EVENT_MONTH<1,prediction_grid$EVENT_MONTH+12,prediction_grid$EVENT_MONTH)
  metrics <- apply(smooth_sample,2,function(x){
    c('ptr'=exp(max(x))/exp(min(x)),
      'peak_month'= month_grid_adjusted[which.max(x)],
      'trough_month'= month_grid_adjusted[which.min(x)])
  })
  ptr_CI <- quantile(metrics['ptr',],probs=c(0.025,0.975))
  month_CI <- apply(metrics[c('peak_month','trough_month'),],1,function(x){
    #detect if vector breaks year cycle (lateest third and first third of the year included)
    if(any(x<4) & any(x>8)){
      x <- ifelse(x>6,x-12,x)
    }
    q_x <- quantile(x,probs=c(0.025,0.975))
    q_x <- ifelse(q_x<1,q_x+12,q_x)
    return(q_x)
  })
  return(list(ptr=ptr_CI,month=month_CI))
}

#' Prepare output for a seasm model object
#'
#' Package the seasonality model to a file to send back the results
#' @param x a named list of seasonality model objects, where the names represent the corresponding the diseases
#' @param ... Not used for this function
#' @param path path to the RData file to be written to disk
#'                    \itemize{
#'                       \item{'-'}{ the data field has been removed, to ensure no individual-level data is included}
#'                       \item{"-"}{ the monthly_counts table has been curated such that only months where 5 or more individuals are diagnosed are included.}
#'                       \item{"-"}{ the monthly counts table contained in both fit and fit_null has been removed}
#'                    }
#'
#'
#' @return No return value, called for side effects. Saves an RData file containing a curated version of the seasm model objects where
#' @seealso \code{\link{seasonality_gam}}, \code{\link{summary.seasm}}, \code{\link{plot.seasm}}
#' @importFrom dplyr tibble bind_rows
#' @export
prepare_output <- function(x,path='output.RData'){
  x_curated <- lapply(x,function(obj){
    obj$data <- NULL
    obj$fit$model <- NULL
    obj$fit_null$model <- NULL
    obj$monthly_counts <- filter(obj$monthly_counts,COUNT>=5)
    return(obj)
  })
  names(x_curated) <- names(x)
  save(x_curated,file=path)
}


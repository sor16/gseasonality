#' Run seasonality model using monthly counts data
#'
#' Runs a GAM which includes an annual and a seasonal term assuming a quasipoisson count distribution
#' @param data data.frame with the diagnosis date for each individual with the disease. It should contain the columns:
#'                    \itemize{
#'                       \item{"ID"}{ a character denoting the ID of the patient}
#'                       \item{"EVENT_DATE"}{ character in the format YYYY-MM-DD denoting the first date the disease is recorded in the registry}
#'                    }
#' @param monthly_counts data.frame with the diagnosis counts per year per month. It should contain the columns:
#'                       \itemize{
#'                          \item{"EVENT_YEAR"}{ Integer denoting the year}
#'                          \item{"EVENT_DATE"}{ Integer denoting the month}
#'                          \item{"COUNT"}{ Integer denoting the number of individuals diagnosed in a particular year and month}
#'                       }
#'                       Either data or monthly_counts must be provided. If both are provided, monthly_counts argument will be ignored and data used to derive the monthly_counts.
#' @param year_start integer denoting the starting year of the analysis.
#' Defaults to NULL in which case the starting year is taken to be the earliest year in data
#' @param year_end integer denoting the final year of the analysis.
#' Defaults to NULL in which case the final year is taken to be the latest year in data
#' @param adjusted boolean denoting whether to perform specific adjustments for months with low hospital usage (July and December). Defaults to FALSE.
#' @param ... Not used for this function
#' @details Details here
#' @return The function returns an object of type "seasm", including the objects:
#'  \item{\code{fit}}{ GAM seasonality model fit object which was inferred using aggregated monthly counts data}
#'  \item{\code{fit_null}}{GAM null model fit object, which does not include a seasonality smooth term, also inferred using aggregated monthly counts data}
#'  \item{\code{fit_table}}{Table with the fitted values for the seasonality and the null model}
#'  \item{\code{seasonality_term}}{Table with the values of the seasonality smooth term from the seasonality model}
#'  \item{\code{annual_term}}{Table with the values of the annual smooth term from the seasonality model}
#'  \item{\code{summary}}{Summary statistics of the seasonality curve}
#'  \item{\code{monthly_counts}}{Diagnosis counts per month, derived from the input data.frame}
#'  \item{\code{data}}{The input data.frame}
#'  \item{\code{year_start}}{input argument or derived from data}
#'  \item{\code{year_end}}{input argument or derived from data}
#'  \item{\code{adjusted}}{input argument}
#' @seealso \code{\link{summary.seasm} and \link{plot.seasm}} for further studying the output
#' @importFrom dplyr mutate filter count summarise group_by
#' @importFrom lubridate ymd year month
#' @importFrom magrittr "%>%"
#' @export
seasonality_gam <- function(data=NULL,monthly_counts=NULL,year_start=NULL,year_end=NULL,adjusted=F){
  if(is.null(data)){
    if(is.null(monthly_counts)){
      stop('Either data or monthly_counts have to be provided')
    }
  }else{
    stopifnot(inherits(data,'data.frame'))
    if(!all(c('ID','EVENT_DATE') %in% names(data))){
      stop('data must include columns ID and EVENT_DATE. Check the help page for more details by writing ?seasonality_gam')
    }
    if(length(unique(data$ID)) != nrow(data)){
      stop('There must be exactly one entry for each ID. Check the help page for more details by writing ?seasonality_gam')
    }
    if(suppressWarnings(any(is.na(ymd(data$EVENT_DATE))))){
      stop('EVENT_DATE column in data must be in the format YYYY-MM-DD. Check the help page for more details by writing ?seasonality_gam')
    }
    data_filtered <- mutate(data,EVENT_DATE=ymd(EVENT_DATE)) %>%
                     mutate(EVENT_YEAR=year(EVENT_DATE),
                            EVENT_MONTH=month(EVENT_DATE),
                            EVENT_DAY=day(EVENT_DATE))
  }
  if(is.null(monthly_counts)){
    if(is.null(data)){
      stop('Either data or monthly_counts have to be provided')
    }
  }else{
    stopifnot(inherits(monthly_counts,'data.frame'))
    if(!all(c('EVENT_YEAR','EVENT_MONTH','COUNT') %in% names(monthly_counts))){
      stop('data must include columns EVENT_YEAR, EVENT_MONTH and COUNT. Check the help page for more details by writing ?seasonality_gam')
    }
    if(!is.null(data)){
      warning('Since data was provided monthly_counts will be ignored. monthly_counts directly generated from data')
    }

  }
  stopifnot(inherits(adjusted,'logical'))
  if(is.null(year_start)){
    year_start <- min(data_filtered$EVENT_YEAR)
  }
  if(is.null(year_end)){
    year_end <- max(data_filtered$EVENT_YEAR)
  }
  check_year_arg(year_start,'year_start')
  check_year_arg(year_end,'year_end')
  year_start <- as.integer(year_start)
  year_end <- as.integer(year_end)
  if(year_end<=year_start){
    stop('year_end must be strictly greater than year_start')
  }
  if(!is.null(data)){
    data_filtered <- filter(data_filtered,EVENT_YEAR>=year_start,EVENT_YEAR<=year_end)
    monthly_counts <- count(data_filtered,EVENT_YEAR,EVENT_MONTH,name = 'COUNT')
  }else{
    monthly_counts <- filter(monthly_counts,EVENT_YEAR>=year_start,EVENT_YEAR<=year_end)
    data_filtered <- NULL
  }
  adjustment <- ifelse(adjusted,'binary','')
  mod_list <- run_seasonality_gam(dat=monthly_counts,adjustment=adjustment,a=1,b=13)
  seasonality_summary <- summarise_seasonality(mod_list=mod_list,monthly_counts=monthly_counts,adjustment=adjustment)
  seasonality_obj <- list()
  attr(seasonality_obj, "class") <- "seasm"
  seasonality_obj$fit <- mod_list$seasonal
  seasonality_obj$fit_null <- mod_list$null
  seasonality_obj$fit_table <- get_fit_table(mod_list=mod_list,year_start=year_start,year_end=year_end)
  seasonality_obj$seasonality_term <- get_seasonality_term(mod=mod_list$seasonal,year_start=year_start,year_end=year_end,adjustment=adjustment)
  seasonality_obj$annual_term <- get_annual_term(mod=mod_list$seasonal,year_start=year_start,year_end=year_end,adjustment=adjustment)
  seasonality_obj$summary <- seasonality_summary
  seasonality_obj$data <- data_filtered
  seasonality_obj$monthly_counts <- monthly_counts
  seasonality_obj$year_start <- year_start
  seasonality_obj$year_end <- year_end
  seasonality_obj$adjusted <- adjusted

  return(seasonality_obj)
}

#' @importFrom stats approx
get_seasonal_spline_adj <- function(seasonal_spline_avg,months){
  approx(x=seasonal_spline_avg$month,
         y=seasonal_spline_avg$avg_seasonal_val,
         xout=months)$y
}

get_month_size <- function(year_start,year_end){
  dates <- seq(ymd(paste0(year_start,'-01-01')),ymd(paste0(year_end,'-12-31')),by='1 day')
  month_size_dat <- tibble(EVENT_YEAR=year(dates),
                           EVENT_MONTH=month(dates),
                           EVENT_DAY=day(dates)) %>%
                    group_by(EVENT_YEAR,EVENT_MONTH) %>%
                    summarise(nr_days=max(EVENT_DAY),.groups = 'drop')
  return(month_size_dat)
}

get_grid <- function(type,year_start,year_end,adjustment='',by_year=0.1,by_month=0.1,seasonal_spline_avg=NULL){
  if(type=='seasonal'){
    grid_dat <-data.frame(EVENT_YEAR=year_start, EVENT_MONTH=seq(0.5,12.5,by=by_month))
  }else if(type=='annual'){
    grid_dat <- data.frame(EVENT_YEAR=seq(year_start,year_end,by=by_year),EVENT_MONTH=1)
  }else{
    stop('Type not recognized')
  }
  #pick arbitrary value for nr_days
  grid_dat$nr_days <- 30
  if(grepl('mean',adjustment)){
    grid_dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,EVENT_MONTH)
  }else if(adjustment=='binary'){
    grid_dat$july <- as.integer(floor(grid_dat$EVENT_MONTH)==7)
    grid_dat$december <- as.integer(floor(grid_dat$EVENT_MONTH)==12 | grid_dat$EVENT_MONTH<1)
  }
  return(grid_dat)
}

#' @importFrom mgcv gam
#' @importFrom lubridate ymd year month day
#' @importFrom dplyr mutate tibble count group_by summarise inner_join
#' @importFrom stats quasipoisson as.formula
#' @importFrom magrittr "%>%"
run_seasonality_gam <- function(dat,a,b,adjustment='',adj_month_length=F,mod_null='adj',seasonal_spline_type='cp',k_seasonal=6,seasonal_spline_avg=NULL){
  month_size_dat <- get_month_size(year_start=min(dat$EVENT_YEAR),year_end=max(dat$EVENT_YEAR))
  dat <- inner_join(dat,month_size_dat,by=c('EVENT_YEAR','EVENT_MONTH'))
  #dat$nr_days <- 30
  k_annual <- 6
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
  month_trough <- ifelse(month_trough<1,month_trough+12,month_trough) # adjust to month range [1,13)
  ptr_est <- val_peak / val_trough
  CI <- get_monte_carlo_PTR_CI(mod=mod_list$seasonal,monthly_counts=monthly_counts,nr_iter=100,adjustment=adjustment,seasonal_spline_avg = seasonal_spline_avg)
  peak_trough <- data.frame(peak=c(month_peak,val_peak),trough=c(month_trough,val_trough))
  idx_cross_over_points <- which(seasonal_component_monthly_est[1:(length(seasonal_component_monthly_est)-1)]*
                                 seasonal_component_monthly_est[2:length(seasonal_component_monthly_est)] < 0)
  cross_over_points <- seasonal_grid$EVENT_MONTH[idx_cross_over_points]
  #note: this can fail if there are two consecutive zeros
  cross_over_list <- list()
  for(i in 1:length(cross_over_points)){
    if(seasonal_component_monthly_est[idx_cross_over_points[i]]>seasonal_component_monthly_est[idx_cross_over_points[i]+1]){
      cross_over_type <- 'high_to_low'
    }else{
      cross_over_type <- 'low_to_high'
    }
    cross_over_list[[cross_over_type]] <- cross_over_points[i]
  }
  rownames(peak_trough) <- c('month','value')
  summary_list <- list(ptr=data.frame(estimate=ptr_est,lower=CI$ptr['2.5%'],upper=CI$ptr['97.5%']),
                       peak_trough=peak_trough,
                       cross_over_points=cross_over_list,
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


check_year_arg <- function(arg,name){
  if(!inherits(arg,'integer')){
    if(inherits(arg,'numeric')){
      if(as.integer(arg)!=arg){
        stop(paste0(name,' must be an integer'))
      }
    }else{
      stop(paste0(name,' must be an integer'))
    }
  }
}


#' @importFrom dplyr as_tibble mutate select
#' @importFrom mgcv gam
get_fit_table <- function(mod_list,year_start,year_end){
  month_size_dat <- get_month_size(year_start=year_start,year_end=year_end)
  dat_grid <- expand.grid(EVENT_MONTH=seq(1,12),EVENT_YEAR=seq(year_start,year_end)) %>%
              as_tibble() %>%
              mutate(july=as.integer(floor(EVENT_MONTH)==7),
                     december=as.integer(floor(EVENT_MONTH)==12 | EVENT_MONTH<1)) %>%
              inner_join(month_size_dat,by=c('EVENT_YEAR','EVENT_MONTH'))
  pred_grid <- mutate(dat_grid,fit_null=predict(mod_list$null,newdata=dat_grid),
                               fit=predict(mod_list$seasonal,newdata=dat_grid))
  return(select(pred_grid,EVENT_YEAR,EVENT_MONTH,fit_null,fit))
}

#' @importFrom dplyr tibble
#' @importFrom mgcv gam
get_seasonality_term <- function(mod,year_start,year_end,adjustment=''){
  months_vec <- seq(0.5,12.5,by=0.01)
  seasonal_pred <- predict(mod,newdata=get_grid(type='seasonal',year_start=year_start,year_end=year_end,adjustment=adjustment,by_month = 0.01),type='terms',unconditional=T,se.fit=T)
  seasonal_pred_dat <- tibble(month=months_vec,
                              est=seasonal_pred$fit[,'s(EVENT_MONTH)'],
                              lower=est - 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'],
                              upper=est + 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'])
  return(seasonal_pred_dat)
}

#' @importFrom dplyr tibble
#' @importFrom mgcv gam
get_annual_term <- function(mod,year_start,year_end,adjustment=''){
  years_vec <- seq(year_start,year_end,by=0.1)
  annual_pred <- predict(mod,newdata=get_grid(type='annual',year_start=year_start,year_end=year_end,adjustment=adjustment),type='terms',unconditional=T,se.fit=T)
  annual_pred_dat <- tibble(year=years_vec,
                            est=annual_pred$fit[,'s(EVENT_YEAR)'],
                            lower=est - 1.96*annual_pred$se.fit[,'s(EVENT_YEAR)'],
                            upper=est + 1.96*annual_pred$se.fit[,'s(EVENT_YEAR)'])
  return(annual_pred_dat)
}

#' Prepare output for a seasm model object
#'
#' Package the seasonality model to a file to send back the results
#' @param x a named list of seasonality model objects, where the names represent the corresponding diseases
#' @param ... Not used for this function
#' @param path path to the RDS file to be written to disk
#'                    \itemize{
#'                       \item{'-'}{ the data field has been removed, to ensure no individual-level data is included}
#'                       \item{"-"}{ the monthly_counts table has been curated such that only months where 5 or more individuals are diagnosed are included.}
#'                       \item{"-"}{ the monthly counts table contained in both fit and fit_null has been removed}
#'                    }
#'
#' @return No return value, called for side effects. Saves an RDS file containing a curated version of the seasm model objects where
#' @seealso \code{\link{seasonality_gam}}, \code{\link{summary.seasm}}, \code{\link{plot.seasm}}
#' @importFrom dplyr mutate tibble bind_rows
#' @export
prepare_output <- function(x,path='output.rds'){
  x_curated <- lapply(x,function(obj){
    obj$data <- NULL
    obj$fit <- NULL
    obj$fit_null <- NULL
    obj$monthly_counts <- mutate(obj$monthly_counts,COUNT=ifelse(COUNT>=5,COUNT,0))
    return(obj)
  })
  names(x_curated) <- names(x)
  saveRDS(x_curated,file=path)
}

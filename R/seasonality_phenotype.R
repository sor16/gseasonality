#' Calculate seasonality phenotype based on seasonality GAM
#'
#' Based on an input model, this function calculates the seasonality phenotype for each individual
#' @param mod object of type "seasm"
#' @param year_start integer denoting the starting year of the analysis
#' @param year_end integer denoting the final year of the analysis
#' @param ... Not used for this function
#' @details Details here
#' @return The function returns a tibble containing two seasonality phenotype:
#'                    \itemize{
#'                       \item{"seasonal_val_binary"}{ Dichotomous value indicating whether an individual was diagnosed during high/low season as determined by mod}
#'                       \item{"seasonal_val_01"}{ The proportional height of the seasonality pattern from mod at the diagnosis date}
#'                    }
#' @seealso \code{\link{seasonality_gam}}
#' @importFrom dplyr mutate filter summarise group_by inner_join select
#' @importFrom lubridate ymd year month day
#' @importFrom mgcv gam
#' @importFrom magrittr "%>%"
#' @importFrom stats approx qnorm
#' @export
get_seasonality_phenotype <- function(mod,year_start,year_end){
  dates <- seq(ymd(paste0(year_start,'-01-01')),ymd(paste0(year_end,'-12-31')),by='1 day')
  num_days_months <- tibble(year=year(dates),month=month(dates),day=day(dates)) %>%
                     group_by(year,month) %>%
                     summarise(num_days=max(day),.groups='drop')
  pheno_dat <- mod$data
  #Convert date to month number
  pheno_dat <- inner_join(pheno_dat,num_days_months,by=c('EVENT_YEAR'='year','EVENT_MONTH'='month')) %>%
               mutate(EVENT_MONTH_DEC=EVENT_MONTH + EVENT_DAY/num_days) %>%
               mutate(EVENT_MONTH_DEC=ifelse(EVENT_MONTH_DEC>12.5,EVENT_MONTH_DEC-12,EVENT_MONTH_DEC))
  month_grid <- get_grid(type='seasonal',year_start=year_start,year_end=year_end)
  seasonal_smooth_term <- predict(mod$fit,newdata=month_grid,type='terms',unconditional=T,se.fit=T)
  pheno_dat <- mutate(pheno_dat,seasonal_val=approx(x=month_grid$EVENT_MONTH,
                                                    y=seasonal_smooth_term$fit[,'s(EVENT_MONTH)'],
                                                    xout=EVENT_MONTH_DEC)$y,
                               seasonal_val_01=(seasonal_val-min(seasonal_val))/(max(seasonal_val)-min(seasonal_val)),
                               seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                               seasonal_val_binary=as.integer(seasonal_val > 0))
  return(select(pheno_dat,ID,EVENT_DATE,EVENT_MONTH_DEC,seasonal_val,seasonal_val_01,seasonal_val_qt,seasonal_val_binary))
}

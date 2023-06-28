#' Calculate seasonality phenotype based on seasonality GAM
#'
#' Based on an input model, this function calculates the seasonality phenotype for each individual
#' @param mod object of type "seasm"
#' @param data Optional data.frame with columns:
#'                    \itemize{
#'                       \item{"ID"}{ a character denoting the ID of the patient}
#'                       \item{"EVENT_DATE"}{ character in the format YYYY-MM-DD denoting the first date the disease is recorded in the registry}
#'                    }
#'              Useful if mod was fitted on an external data set.
#' @param year_start integer denoting the starting year of the analysis
#' @param year_end integer denoting the final year of the analysis
#' @param adjusted boolean denoting whether to perform specific adjustments for months with low hospital usage (July and December). Defaults to FALSE.
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
get_seasonality_phenotype <- function(mod,data=NULL,year_start=NULL,year_end=NULL,adjusted=F){
  stopifnot(inherits(mod,'seasm'))
  if(!is.null(data)){
    stopifnot(inherits(data,'data.frame'))
    if(!all(c('ID','EVENT_DATE') %in% names(data))){
      stop('data must include columns ID and EVENT_DATE. Check the help page for more details by writing ?get_seasonality_phenotype')
    }
    if(length(unique(data$ID)) != nrow(data)){
      stop('There must be exactly one entry for each ID. Check the help page for more details by writing ?get_seasonality_phenotype')
    }
    if(suppressWarnings(any(is.na(ymd(data$EVENT_DATE))))){
      stop('EVENT_DATE column in data must be in the format YYYY-MM-DD. Check the help page for more details by writing ?get_seasonality_phenotype')
    }
    if(!is.null(mod$data)){
      warning('Since data is provided, it will be used instead of the data in the seasonality model object.')
    }
    pheno_dat <- mutate(data,EVENT_DATE=ymd(EVENT_DATE)) %>%
                 mutate(EVENT_YEAR=year(EVENT_DATE),
                        EVENT_MONTH=month(EVENT_DATE),
                        EVENT_DAY=day(EVENT_DATE))
  }else{
    if(is.null(mod$data)){
      stop('No individual-level data has been provided. Can either be included in mod or provided using the data argument. For more information, write ?seasonality_phenotype')
    }
    pheno_dat <- mod$data
  }
  if(is.null(year_start)){
    year_start <- mod$year_start
  }
  if(is.null(year_end)){
    year_end <- mod$year_end
  }
  check_year_arg(year_start,'year_start')
  check_year_arg(year_end,'year_end')
  year_start <- as.integer(year_start)
  year_end <- as.integer(year_end)
  if(year_end<=year_start){
    stop('year_end must be strictly greater than year_start')
  }
  month_size_dat <- get_month_size(year_start,year_end)
  #Convert date to month number
  pheno_dat <- inner_join(pheno_dat,month_size_dat,by=c('EVENT_YEAR','EVENT_MONTH')) %>%
               mutate(EVENT_MONTH_DEC=EVENT_MONTH + EVENT_DAY/nr_days) %>%
               mutate(EVENT_MONTH_DEC=ifelse(EVENT_MONTH_DEC>12.5,EVENT_MONTH_DEC-12,EVENT_MONTH_DEC))
  #Create phenotypes based on seasonal smooth term
  pheno_dat <- mutate(pheno_dat,seasonal_val=approx(x=mod$seasonality_term$month,
                                                    y=mod$seasonality_term$est,
                                                    xout=EVENT_MONTH_DEC)$y,
                               seasonal_val_01=(seasonal_val-min(seasonal_val))/(max(seasonal_val)-min(seasonal_val)),
                               seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                               seasonal_val_binary=as.integer(seasonal_val > 0))
  return(select(pheno_dat,ID,EVENT_DATE,EVENT_MONTH_DEC,pheno_qt=seasonal_val_qt,pheno_binary=seasonal_val_binary))
}

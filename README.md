
<!-- README.md is generated from README.Rmd. Please edit that file -->

## gseasonality - Genetics of seasonality <img src="man/figures/logo.png" align="right" alt="" width="140" />

This software package enables assessing the role of genetics in seasonal
disease risk. In other words, assessing whether genetic variants are
associated with a tendency to be diagnosed with a disease at a certain
time of year. There are two main functions in the package:

`seasonality_gam` - Infers the seasonality model using disease diagnosis
dates

`get_seasonality_phenotype` - Assigns a seasonality phenotype to each
individual based on the seasonality model.

## Installation

If you are connected to internet, the package can be installed by
writing

``` r
devtools::install_github("sor16/gseasonality")
```

For an environment not connected to internet, download a zip file of the
package, ship it to your environment and write in the RStudio console:

``` r
devtools::install_local('path/to/package/zip/file',type='source',repos=NULL)
```

## Getting started

The gseasonality package is easy to use. The first step is to fit the
seasonality model. It comes with a simulated data set `diagnosis_data`
which demonstrates the desired format of the data for the seasonality
model. It has to have columns with names ID and EVENT_DATE. To run the
seasonality model, write

``` r
library(gseasonality)
seasonality_mod <- seasonality_gam(diagnosis_data)
```

The output of the function is type ‘seasm’ and `summary` and `plot`
functions can be applied to it to explore and visualize the results,
e.g.:

``` r
summary(seasonality_mod)
#> 
#> PTR:
#>  estimate lower-2.5% upper-97.5%
#>      1.39       1.32        1.48
#> 
#> Characteristics:
#>       peak trough
#> month 9.09   3.02
#> value 1.16   0.83
#> 
#> Deviance explained: 37.28%
#> Dispersion: 0.97
```

and

``` r
plot(seasonality_mod)
```

![](man/figures/unnamed-chunk-6-1.png)<!-- --> (plot not shown).

Next, to create the seasonality phenotypes based on the seasonality
model fit, we write:

``` r
seasonality_pheno <- get_seasonality_phenotype(seasonality_mod)
head(seasonality_pheno)
#>   ID EVENT_DATE EVENT_MONTH_DEC seasonal_val_qt seasonal_val_binary
#> 1  1 2019-01-31        2.000000      -1.0431363                   0
#> 2  2 2005-07-07        7.225806       0.4460416                   1
#> 3  3 2009-07-14        7.451613       0.5503623                   1
#> 4  4 2011-03-07        3.225806      -1.9565748                   0
#> 5  5 2009-12-05       12.161290      -0.2005404                   1
#> 6  6 2012-06-01        6.033333      -0.1315310                   1
```

The output is a tibble data frame containing:

- `seasonal_val_binary` - Binary seasonality phenotype. 1 denotes that
  the individual was diagnosed in the season of high incidence and 0
  denotes that the individual was diagnosed in the season of low
  incidence

- `seasonal_val_qt` - Continuous seasonality phenotype value. This
  denotes the extent of the seasonality curve when an individual is
  diagnosed with a disease. The phenotype has been quantile normalized.

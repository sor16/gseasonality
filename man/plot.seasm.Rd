% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasm_methods.R
\name{plot.seasm}
\alias{plot.seasm}
\title{Plot method for seasonality GAM on monthly counts}
\usage{
\method{plot}{seasm}(x, ..., type = "panel")
}
\arguments{
\item{...}{other plotting parameters (not used in this function)}

\item{type}{a character denoting what type of plot should be drawn. Defaults to "rating_curve". Possible types are
\itemize{
   \item{"fit"}{ to plot the GAM fit to the monthly counts}
   \item{"seasonality"}{ to plot the seasonal smooth term in th GAM model}
   \item{"annual"}{ to plot the annual smooth term in th GAM model}
   \item{"panel"}{ to plot all three plots above (default)}
}}

\item{object}{an object of class "seasm".}
}
\value{
No return value, called for side effects.
}
\description{
Visualize different aspects of the seasonality GAM
}
\seealso{
\code{\link{seasonality_gam}}, \code{\link{summary.seasm}}, \code{\link{autoplot.seasm}}
}

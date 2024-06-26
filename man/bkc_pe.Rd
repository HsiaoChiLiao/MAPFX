% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/022.bkc_pe.R
\name{bkc_pe}
\alias{bkc_pe}
\title{Background correction for the well-specific markers (PE)}
\usage{
bkc_pe(paths, pe.lower.quantile = 0.1, pe.min.quantile = 0.01, plots = TRUE)
}
\arguments{
\item{paths}{a vector of characters of paths to store intput, intermediary results, outputs...}

\item{pe.lower.quantile}{the cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise for infinity markers.}

\item{pe.min.quantile}{the cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers for infinity markers.}

\item{plots}{logical; if TRUE (default), produce scatter plots for pre- and post- background adjusted infinity markers (calibrated values on y-axis and raw values on x-axis).}
}
\value{
Background noise corrected infinity markers and graphs if specified
}
\description{
This function has been designed to do background correction for the well-specific markers (PE) by using normal-exponential convolution model.
}
\details{
Generating the calibrated measurements and save to bkc.pe_mt.rds file, and visualising the result with the scatter plots in the output directory.
}
\author{
Hsiao-Chi Liao and Agus Salim
}

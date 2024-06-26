% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/040.mkImputeMT.R
\name{mkImputeMT}
\alias{mkImputeMT}
\title{Making the input for imputation}
\usage{
mkImputeMT(paths)
}
\arguments{
\item{paths}{a vector of characters of paths to store intput, intermediary results, outputs...}
}
\value{
Combined normalised backbone and infinity markers for imputation
}
\description{
This function has been designed to combine the normalised backbone measurements and the normalised PE markers for later imputation.
}
\details{
Generating the combined data and saving to impu.input_log.mt.rds (on log scale) in the output directory.
}
\author{
Hsiao-Chi Liao
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_for_bma.r
\name{convert_for_bma}
\alias{convert_for_bma}
\title{Convert LOESS Summary List to BMA-Compatible Data Frame}
\usage{
convert_for_bma(loess_summaries)
}
\arguments{
\item{loess_summaries}{A list object returned by \code{smooth_occti_outputs()},
which must contain a \code{loess_summaries} element (a named list of data frames per species).}
}
\value{
A data frame with columns:
\describe{
\item{species}{Species name from the list names}
\item{year}{Year of the estimate}
\item{index}{Smoothed occupancy index (mean)}
\item{se}{Standard error of the smoothed index}
}
}
\description{
This function converts the LOESS summary outputs from \code{smooth_occti_outputs()} into a format
suitable for Bayesian Model Averaging (BMA) input. It binds all species data into a single
data frame and renames key columns to match expected BMA input structure.
}

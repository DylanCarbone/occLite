% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_occ_ms.R
\name{fit_occ_ms}
\alias{fit_occ_ms}
\title{Wrapper function to fit occupancy models to multiple species.
#' Function to fit the occupancy models}
\usage{
fit_occ_ms(
  splist,
  obdata,
  occformula = "~North+I(North^2)+East+I(East^2)",
  detformula = "~logLL+SEAS",
  covnames = c("East", "North"),
  seed = NULL,
  minyear = NULL,
  maxyear = NULL,
  trendyears = NULL,
  nstart = 3,
  nrun = 1,
  allsites = NULL,
  qval = 0.025,
  prev_start = NULL,
  outputdir = NULL,
  printprogress = FALSE,
  engine = "C"
)
}
\arguments{
\item{splist}{Vector of target species to model}

\item{obdata}{Data frame containing species occurrence records with the following columns: species, Date, gridref, Year, Week, Month, and optionally covnames (see below) and listL}

\item{occformula}{Formula for occupancy probability}

\item{detformula}{Formula for detection probability}

\item{covnames}{Vector of covariate names in obdata}

\item{seed}{Option random seed value to set}

\item{minyear}{Optional filter to first year of interest}

\item{maxyear}{Optional filter to last year of interest}

\item{trendyears}{Vector of start years for trend estimation. If \code{trendyears = NULL} then no trends will be calculated.}

\item{nstart}{Number of starting values to run. Default \code{nstart = 3}.}

\item{nrun}{Optional number of runs (since starting values are based on preceding year). Default \code{nrun = NULL} produces one run.}

\item{allsites}{Optional data frame of sites for which the occupancy index will be calculated for.}

\item{qval}{Quantile value to filter records to months where the species was observed. Default \code{qval = 0.025}.}

\item{prev_start}{Provide starting values e.g. based on outputs of a previous run.}

\item{outputdir}{Optional directory to save output files to.}

\item{printprogress}{print the progress of the run (only available for non-parallel option)}

\item{engine}{Choose the engine used by unmarked.}
}
\value{
A list containing various outputs
}
\description{
Wrapper function to fit occupancy models to multiple species.
#' Function to fit the occupancy models
}

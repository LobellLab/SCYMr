% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SentinelCalibration.R
\name{calibrateData}
\alias{calibrateData}
\title{Apply SR calibration to Sentinel TOA data}
\usage{
calibrateData(sentineldf, toaColName, int, slope)
}
\arguments{
\item{sentineldf}{Data.frame containing paired Sentinel TOA and SR values for a band or index}

\item{toaColName}{column name for Sentinel TOA values}

\item{slope}{slope from linear calibration model}

\item{intercept}{intercept from linear calibration model}
}
\description{
Applies the calibration to sentinel data based on a single intercept and slope
}
\examples{
# apply calibration
calibrated <- calibrateData(sentinelData_df, 'Sentinel_TOA', 0.05, 1.903)
}
\keyword{Calibration}
\keyword{Sentinel}

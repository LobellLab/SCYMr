% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scym_core.R
\name{runScymCore_peak}
\alias{runScymCore_peak}
\title{Run SCYM Peak: basic 1-window SCYM function based on US maize variables}
\usage{
runScymCore_peak(doy, vi1, p_hinge, v_hinge, Augtmax, JJA_radn,
  apsimtable = ctable)
}
\arguments{
\item{doy}{column name for character vector denoting observation day of year}

\item{vi1}{column name for vi observation at peak}

\item{p_hinge}{column name for precip hinge value}

\item{v_hinge}{column name for vpd hinge value}

\item{Augtmax}{column name for augmaxt}

\item{apsimtable}{data.frame with coefficient table to use}
}
\description{
Function intended to be applied by row to a df of points + attributes
}
\examples{
# to apply to data frame using purrr (better than rowwise, etc)
ls2017_yield <- ls2017 \%>\%
  mutate(biomass_scymr = pmap_dbl(list(doy = DOY, vi1 = GCVI_max,
                                       p_hinge=phinge, v_hinge=vhinge,
                                       Augtmax = Augmaxt, JJA_radn = JJAradn),
                                       .f = runScymCore, apsimtable = ctable),
         yield_SCYMr_tha = biomass_scymr * 0.45)
}
\keyword{SCYM}
\keyword{US}
\keyword{core,}
\keyword{default}
\keyword{maize}

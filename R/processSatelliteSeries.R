#' Extract the maximum index value by crop window for each point from raw csv files
#'
#' Reduces a time series of satellite observations for each point to a maximum
#' value for the specified windows. Requires readr, lubridate, and dplyr.
#' Initial setup is for
#' exports from Jill's GEE script '02.1x_satelliteSampler_forScymOffline'
#' @param df_file Data.frame containing a veg index of interest
#' @param datasource "Landsat" or "Sentinel"
#' @param maskClouds should cloud qa be applied? defaults to FALSE
#' @param w1_doy1 first day of first window; default is US maize
#' @param w1_doy2 last day of first window; default is US maize
#' @param w2_doy1 first day of second window; default is US maize
#' @param w2_doy2 last day of second window; default is US maize
#' @keywords SCYM preparation, GCVI; default is US maize
#' @export
#' @examples
#' # to apply to a list of files:
#' landsatFiles <- list.files(gDrive, pattern=landsatBaseName, full.names = TRUE)
#' # check file sizes and drop empty exports
#' landsatFilesb <- landsatFiles[sapply(landsatFiles, file.size) > 2]
#' # load/get window GCVIs - default windows
#' landsat_gcvi <- purrr::map_df(landsatFilesb, getWindowGCVIs_fromFile,
#'                                datasource = "Landsat", maskClouds = TRUE)

getWindowGCVIs_fromFile <- function(df_file, datasource, maskClouds = FALSE,
                                    w1_doy1 = 165, w1_doy2 = 200, w2_doy1 = 201,
                                    w2_doy2 = 240){
  library(readr)
  library(dplyr)
  library(lubridate)

  # load and drop points outside windows
  griddf <- readr::read_csv(df_file, guess_max = 1000) %>%
    dplyr::select(-c(`system:index`, `.geo`, yield_scym, yield_tha, state2)) %>%
    mutate(date = ymd(date),
           doy = yday(date)) %>%
    # assign window
    mutate(window = case_when(doy < w1_doy1 ~ 0,
                              doy >= w1_doy1 & doy <= w1_doy2 ~1,
                              doy >= w2_doy1 & doy <= w2_doy2 ~2,
                              doy > w2_doy2 ~ 0)) %>%
    filter(window != 0)

  # mask clouds if flagged
  if(maskClouds){
    if(datasource == 'Landsat'){
      griddf <- griddf %>%
        filter(pxqa_clear == 1)
    }
    if(datasource == 'Sentinel'){
      griddf <- griddf %>%
        filter(QA60_DECODED == 1)
    }
  }

  # get max gcvi per window per point
  maxes <- griddf %>%
    mutate(GCVI = (NIR/GREEN)-1) %>%
    group_by(pointID, window, year) %>%
    slice(which.max(GCVI)) %>%
    dplyr::select(c(pointID, year, fips, state, gridID, granularID, window, doy, GCVI))

  # yes, this is oddly inefficient
  maxes_wider <- maxes %>%
    tidyr::spread(., key = window, value = doy) %>%
    rename(doy1 = `1`,
           doy2 = `2`)
  maxwin1 <- maxes_wider %>%
    dplyr::select(-c(doy2)) %>%
    rename(gcvi1 = GCVI) %>%
    tidyr::drop_na()
  maxwin2 <- maxes_wider %>%
    dplyr::select(c(pointID, year, doy2, GCVI)) %>%
    rename(gcvi2 = GCVI) %>%
    tidyr::drop_na()
  # combine
  max_wide <- maxwin1 %>%
    left_join(maxwin2) %>%
    mutate(Dates = paste0(doy1, '_', doy2))

  return(max_wide %>% ungroup())
}


#' Extract the maximum index value and DOY for each point from raw csv files
#'
#' Reduces a time series of satellite observations for each point to a maximum
#' value for the specified window. Requires readr, lubridate, and dplyr.
#' Default window is June - August.
#' Initial setup is for
#' exports from Jill's GEE script '02.1x_satelliteSampler_forScymOffline'
#' @param df_file Data.frame containing a veg index of interest
#' @param datasource "Landsat" or "Sentinel"
#' @param maskClouds should cloud qa be applied? defaults to FALSE
#' @param w1_doy1 first day of first window; default is June 1
#' @param w1_doy2 last day of first window; default is Sept 1
#' @keywords SCYM preparation, peak vi
#' @export
#' @examples
#' # to apply to a list of files:
#' landsatFiles <- list.files(gDrive, pattern=landsatBaseName, full.names = TRUE)
#' # check file sizes and drop empty exports
#' landsatFilesb <- landsatFiles[sapply(landsatFiles, file.size) > 2]
#' # load/get window GCVIs - default windows
#' landsat_gcvi <- purrr::map_df(landsatFilesb, getPeakGCVIs_fromFile,
#'                                datasource = "Landsat", maskClouds = TRUE)

getPeakGCVIs_fromFile <- function(df_file, datasource, maskClouds = FALSE,
                                    w1_doy1 = 152, w1_doy2 = 245){
  library(readr)
  library(dplyr)
  library(lubridate)

  # load and drop points outside windows
  griddf <- readr::read_csv(df_file, guess_max = 1000) %>%
    dplyr::select(-c(`system:index`, `.geo`, yield_scym, yield_tha, state2)) %>%
    mutate(date = ymd(date),
           doy = yday(date)) %>%
    # filter for window
    filter(doy >= w1_doy1 & doy <= w1_doy2)

  # mask clouds if flagged
  if(maskClouds){
    if(datasource == 'Landsat'){
      griddf <- griddf %>%
        filter(pxqa_clear == 1)
    }
    if(datasource == 'Sentinel'){
      griddf <- griddf %>%
        filter(QA60_DECODED == 1)
    }
  }

  # get max gcvi per window per point
  maxes <- griddf %>%
    mutate(GCVI = (NIR/GREEN)-1) %>%
    group_by(pointID, year) %>%
    slice(which.max(GCVI)) %>%
    dplyr::select(c(pointID, year, fips, state, doy, GCVI))

  return(maxes)
}

#' Fit harmonics recursivelys to satellite time series
#'
#' Runs a 2-term harmonic regression for the specified number of iterations.
#' Written by Calum You with modification by Jake Campolo. Output includes
#' harmonic values at original data observations, but not coefficients.
#' @param ts Data.frame containing columns date, dependent variable, and any
#' grouping variables.
#' @param iterations Number of recursive harmonics to fit
#' @param origin reference time; format should match df date column
#' @param dependent column name for dependent variable
#' @param dt_col unquoted column name for date column. Default is DATE
#' @param omega default is 1
#' @param ... additional grouping variables; unquoted
#' @keywords harmonics
#' @export
#' @examples
#' # apply harmonics by pointID
#' library(dplyr)
#' library(lubridate)
#' library(stringr)
#' library(broom)
#'
#' df_fitted <- df_satelliteSeries %>%
#'   addHarmonicFits(
#'     iterations = 4,
#'     origin = ymd("2018-01-01"),
#'     dependent = "GCVI",
#'     pointID
#'  )

addHarmonicFits <- function(ts, iterations, origin, dependent, ..., dt_col = DATE, omega = 1){
  library(dplyr)
  library(stringr)
  library(broom)
  library(lubridate)
  dt_col <- enquo(dt_col)
  groupers <- enquos(...)
  dependent = ensym(dependent)

  hdt_ts <- ts %>%
    add_hdt(origin, dt_col = !!dt_col, omega = omega) %>%
    mutate(raw = !!dependent)

  for (i in 1:iterations){
    this_fit = rlang::sym(str_c("Fit", i))
    hdt_ts <- hdt_ts %>%
      add_augment(dependent = "raw", !!!groupers) %>% # Get fitted values
      mutate(
        raw = pmax(fitted, raw), # Make new raw values for next iteration
        !!this_fit := fitted # Name this iterations fitted values
      ) %>%
      select(!!!groupers, raw, matches("Fit\\d"), !!dt_col, t, p1, p2, p3, p4)
  }
  return(hdt_ts)
}

# helpers for addHarmonicFits
add_hdt <- function(ts, origin, dt_col = DATE, omega = 1){

  dt_col <- enquo(dt_col)
  ts %>% mutate(
    daynum = case_when(leap_year(!!dt_col) ~ 366,
                       !leap_year(!!dt_col) ~ 365),
    t = as.numeric(!!dt_col - origin) / daynum,
    p1 = sin(2 * pi * t * omega),
    p2 = cos(2 * pi * t * omega),
    p3 = sin(4 * pi * t * omega),
    p4 = cos(4 * pi * t * omega)
  )
}

add_models <- function(ts, dependent, ...){
  groupers <- enquos(...)
  equation <- dependent %>%
    str_c("~ p1 + p2 + p3 + p4") %>%
    as.formula()
  ts %>%
    group_by(!!!groupers) %>%
    do(model = lm(equation, data = .))
}

add_augment <- function(ts, dependent, ..., dt_col = DATE){
  dt_col <- enquo(dt_col)
  groupers <- enquos(...)
  dep_sym <- ensym(dependent)

  models <- ts %>%
    add_models(dep_sym, ...) %>%
    augment(model)

  grp_str <- groupers %>% toString %>% str_remove_all("~") %>% str_split(",") %>% flatten_chr

  ts %>%
    left_join(models, by = c(grp_str, dependent, "p1", "p2", "p3", "p4")) %>%
    select(!!!groupers, !!dep_sym, fitted = .fitted, !!dt_col, t, p1, p2, p3, p4, matches("Fit\\d"))
}


#' Get final recursive harmonics model coefficients
#'
#' Runs a 2-term harmonic regression for the specified number of iterations.
#' Modifies getHarmonicFits to return the coefficients instead of fitted values
#' @param ts Data.frame containing columns date, dependent variable, and any
#' grouping variables.
#' @param iterations Number of recursive harmonics to fit
#' @param origin reference time; format should match df date column
#' @param dependent column name for dependent variable
#' @param dt_col unquoted column name for date column. Default is DATE
#' @param omega default is 1
#' @param ... additional grouping variables; unquoted
#' @keywords harmonics
#' @export
#' @examples
#' # apply harmonics by pointID
#' library(dplyr)
#' library(lubridate)
#' library(stringr)
#' library(broom)
#'
#' coefficients <- df_satelliteSeries %>%
#'   getRecursiveCoeffs(iterations = numIterations,
#'                      origin = ymd("2018-01-01"),
#'                      dependent = "GCVI",
#'                       omega = 1.5, pointID)

getRecursiveCoeffs <- function(ts, iterations, origin, dependent0, ..., dt_col = DATE, omega = 1){
  library(dplyr)
  library(stringr)
  library(broom)
  dt_col <- enquo(dt_col)
  groupers <- enquos(...)
  dependent = ensym(dependent0)

  hdt_ts <- ts %>%
    add_hdt(origin, dt_col = !!dt_col, omega = omega) %>%
    mutate(raw = !!dependent)

  # get the fitted values
  for (i in 1:iterations){
    this_fit = rlang::sym(str_c("Fit", i))
    hdt_ts <- hdt_ts %>%
      add_augment(dependent = "raw", !!!groupers) %>% # Get fitted values
      mutate(
        raw = pmax(fitted, raw), # Make new raw values for next iteration
        !!this_fit := fitted # Name this iterations fitted values
      ) %>%
      select(!!!groupers, raw, matches("Fit\\d"), !!dt_col, t, p1, p2, p3, p4)
  }

  # rerun last harmonic to get coeffs
  coeffRenamer <- data.frame(term = c('(Intercept)','p1','p2','p3','p4'),
                             independents = c('constant','sin','cos','sin2','cos2'),
                             stringsAsFactors = FALSE)

  lastFit <- paste0('Fit',iterations)
  equationLast <- paste0(lastFit, ' ~ p1 + p2 + p3 + p4')

  bestModels <- hdt_ts %>%
    group_by(!!!groupers) %>%
    do(model = lm(equationLast, data = .))

  # get coefficients
  coefficients <- bestModels %>%
    tidy(model) %>%
    left_join(coeffRenamer) %>%
    dplyr::select(-c(std.error, statistic, p.value, term)) %>%
    tidyr::spread(., key = independents, value = estimate)

  return(coefficients)
}


#' Predict value from harmonics
#'
#' @param doy doy for prediction
#' @param constant
#' @param a1
#' @param a2
#' @param b1
#' @param b2
#' @param omega default is 1.5
#' @keywords harmonics, prediction
#' @export
#' @examples

predictHarmonics <- function(doy, constant, a1, a2, b1, b2, omega = 1.5){
  t <- doy/365
  prediction <- constant + a1*cos(2*t*pi*omega) + b1*sin(2*t*pi*omega) +
    a2*cos(4*t*pi*omega) + b2*sin(4*t*pi*omega)
}



#the red edge chlorophyl index comes from various Gitelson papers (originally in 2003!?)
#pass in VI option: GCVI, CIr, NDIr, NIRv
#creates biweekly time series of the VI, taking max observation in each 2 week window
getBiweeklyData_fromFile <- function(df_file, year, VI, TOAcalib = TRUE) {
  library(dplyr)
  library(tidyverse)
  library(reshape2)
  #create a sequence of two week intervals
  date_intervals <- seq(ymd(paste0(year, '-04-01', sep="")),ymd(paste0(year, '-10-30', sep = "")), by = '2 week')
  #read in file
  griddf <- readr::read_csv(df_file, guess_max = 1000) %>%
    dplyr::select(-c(`system:index`, `.geo`, yield_scym, yield_tha, state2)) %>%
    mutate(date = ymd(date),
           doy = yday(date))
  if (VI == 'GCVI'){
    griddf <- griddf %>% mutate(GCVI = (NIR/GREEN)-1)
  } else if (VI == 'NIRv') {
    griddf <- griddf %>% mutate(NIRv = NIR * (NIR - RED)/(NIR + RED))
  } else if (VI == "CIr") {
    griddf <- griddf %>% mutate(CIr =  (NIR/RDED2) - 1)
  } else if (VI == "NDIr") {
    griddf <- griddf %>% mutate(NDIr =  (RDED1 - RED)/(RDED1 + RED))
  }

  #for each 2 week interval, slice by max GCVI
  hold_maxes <- list()
  for (i in 1:(length(date_intervals)-1)) {
    #set up interval
    interval <- lubridate::interval(date_intervals[i], date_intervals[i+1])
    if (VI == "GCVI") {
      maxObs <- griddf %>% filter(date %within% interval) %>% group_by(pointID) %>%
        slice(which.max(GCVI))
    } else if (VI == "NIRv") {
      maxObs <- griddf %>% filter(date %within% interval) %>% group_by(pointID) %>%
        slice(which.max(NIRv))
    } else if (VI == "CIr") {
      maxObs <- griddf %>% filter(date %within% interval) %>% group_by(pointID) %>%
        slice(which.max(NIR))
    } else if (VI == "NDIr") {
      maxObs <- griddf %>% filter(date %within% interval) %>% group_by(pointID) %>%
        slice(which.max(RDED1))
    }
    hold_maxes[[i]] <- maxObs %>% mutate(period = i)
  }
  maxVI_biweekly <- do.call('rbind', hold_maxes)
  #adjust GCVI
  if (VI == "GCVI") {
    maxVI_biweekly <- calibrateData(maxVI_biweekly, 'GCVI', int = -0.05, slope = 1.903) %>% select(-GCVI) %>%
      rename("GCVI" = "GCVI_adj")
    maxVI_biweekly_wide <- maxVI_biweekly %>% select(c(pointID, period, fips, granularID, state, year, GCVI))
  } else if (VI == "NIRv") {
    maxVI_biweekly <- calibrateData(maxVI_biweekly, 'NIRv', int = 63.21, slope = 1.215) %>% select(-NIRv) %>%
      rename("NIRv" = "NIRv_adj")
    maxVI_biweekly_wide <- maxVI_biweekly %>% select(c(pointID, period, fips, granularID, state, year, NIRv))
  } else if (VI == "CIr") {
     maxVI_biweekly <- calibrateData(maxVI_biweekly, 'CIr', int = 0.05817, slope = 0.98996) %>% select(-CIr) %>%
     rename("CIr" = "CIr_adj")
    maxVI_biweekly_wide <- maxVI_biweekly %>% select(c(pointID, period, fips, granularID, state, year, CIr))
  } else if (VI == "NDIr") {
     maxVI_biweekly <- calibrateData(maxVI_biweekly, 'NDIr', int = 0.01327, slope = 1.62311) %>% select(-NDIr) %>%
     rename("NDIr" = "NDIr_adj")
    maxVI_biweekly_wide <- maxVI_biweekly %>% select(c(pointID, period, fips, granularID, state, year, NDIr))
  }
  #melt and cast
  maxVI_biweekly_wide <- dcast(melt(maxVI_biweekly_wide, id.vars = c("pointID", "period", "fips", "granularID", "state", "year")),
                               pointID + fips + granularID + state + year ~ variable + period)
  return(maxVI_biweekly_wide)
}


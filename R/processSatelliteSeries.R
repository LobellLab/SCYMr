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
    dplyr::select(c(pointID, year, fips, state, window, doy, GCVI))

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

  return(max_wide)
}


#' Fit harmonics recursivelys to satellite time series
#'
#' Runs a 2-term harmonic regression for the specified number of iterations
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
    t = as.numeric(!!dt_col - origin) / 365,
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




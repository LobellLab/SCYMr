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


#' Run SCYM core: basic 2-window SCYM function based on US maize variables
#'
#' Function intended to be applied by row to a df of points + attributes
#' @param days column name for character vector denoting observation days for both windows in format '180_210'
#' @param vi1 column name for vi observation in window 1
#' @param vi2 column name for vi observation in window 2
#' @param p_hinge column name for precip hinge value
#' @param v_hinge column name for vpd hinge value
#' @param Augtmax column name for augmaxt
#' @param apsimtable data.frame with coefficient table to use
#' @keywords SCYM core, default US maize
#' @export
#' @examples
#' # to apply to data frame using purrr (better than rowwise, etc)
#' ls2017_yield <- ls2017 %>%
#'   mutate(biomass_scymr = pmap_dbl(list(days = Dates, vi1 = gcvi1, vi2 = gcvi2,
#'                                        p_hinge=phinge, v_hinge=vhinge,
#'                                        Augtmax = Augmaxt, JJA_radn = JJAradn),
#'                                        .f = runScymCore, apsimtable = ctable),
#'          yield_SCYMr_tha = biomass_scymr * 0.45)


runScymCore <- function(days, vi1, vi2, p_hinge, v_hinge, Augtmax,
                        JJA_radn, apsimtable = ctable){
  # get coeffs
  coeffs <- apsimtable[apsimtable$Dates == days,]

  # yield
  biomass <- coeffs$Intcept + (vi1 * coeffs$gcvi1) + (vi2 *coeffs$gcvi2) +
    (p_hinge * coeffs$phinge) + (v_hinge*coeffs$vhinge) +
    (Augtmax * coeffs$Augmaxt) + (JJA_radn * coeffs$JJAradn)

  return(biomass)
}




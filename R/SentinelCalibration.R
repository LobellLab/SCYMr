#' Get Sentinel TOA to Sentinel SR Calibration
#'
#' Calculates a linear equation between Sentinel SR and TOA. Returns the
#' intercept and slope as a data frame, along with the r^2
#' @param sentineldf Data.frame containing paired Sentinel TOA and SR values for a band or index
#' @param toaColName column name for Sentinel TOA values
#' @param srColName column name for Senintel SR values
#' @keywords Sentinel Calibration
#' @export
#' @examples
#' # get linear coefficients
#' coefficients <- getCalibrationCoeffs(sentinelData_df, 'Sentinel_TOA', 'Sentinel_SR')

getCalibrationCoeffs <- function(sentineldf, toaColName, srColName){
  # run lm and extract coefficients
  formula <- paste0(srColName, '~', toaColName)
  s_lm <- lm(formula, data = sentineldf)
  lm.int <- s_lm$coefficients[[1]]
  lm.slope <- s_lm$coefficients[[2]]

  outdf <- data.frame(intercept = lm.int,
                      slope = lm.slope,
                      r2 = summary(s_lm)$r.squared)
  return(outdf)
}

#' Apply SR calibration to Sentinel TOA data
#'
#' Applies the calibration to sentinel data based on a single intercept and slope
#' @param sentineldf Data.frame containing paired Sentinel TOA and SR values for a band or index
#' @param toaColName column name for Sentinel TOA values
#' @param intercept intercept from linear calibration model
#' @param slope slope from linear calibration model
#' @keywords Sentinel Calibration
#' @export
#' @examples
#' # apply calibration
#' calibrated <- calibrateData(sentinelData_df, 'Sentinel_TOA', 0.05, 1.903)

calibrateData <- function(sentineldf, toaColName, int, slope){
  outColName <- paste0(toaColName, '_adj')
  sentineldf[,outColName] <- int + (slope * sentineldf[,toaColName])
  return(sentineldf)
}

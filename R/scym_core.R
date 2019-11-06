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




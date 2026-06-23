#' Determine linear range status for measured samples
#'
#' @description
#' Assigns each biological sample a linear-range status based on calibration
#' curve boundaries derived from dilution/concentration experiments.
#'
#' For each sample intensity (`y`), the function checks whether it lies:
#' - within the calibrated linear range
#' - below the lower limit of linearity (BLOL)
#' - above the upper limit of linearity (ULOL)
#'
#' Additionally, samples failing the R² quality criterion (`aboveR2 == FALSE`)
#' are flagged as not reliable for linear range interpretation.
#'
#' @param dats Data table containing biological sample information.
#' Must include at least:
#' \itemize{
#'   \item IDintern (unique sample identifier)
#'   \item groupIndices (link to calibration curve)
#'   \item y (measured response variable)
#' }
#'
#' @param datCal Data table containing calibration curve metadata.
#' Must include:
#' \itemize{
#'   \item groupIndices (matching key to `dats`)
#'   \item RangeStartY (lower linear range boundary)
#'   \item RangeEndY (upper linear range boundary)
#'   \item aboveR2 (quality flag for regression model)
#'   \item LRFlag (optional precomputed linear-range flag)
#' }
#'
#' @param y Character. Name of response variable column in `dats`
#' (typically measured intensity or signal).
#'
#' @return
#' Data table identical to `dats` with additional columns:
#' \itemize{
#'   \item Status_LR: classification of linear range status
#'     ("BLOL", "ULOL", TRUE/FALSE)
#'   \item LRFlag: inherited calibration flag
#'   \item RangeStartY: lower boundary of linear range
#'   \item RangeEndY: upper boundary of linear range
#' }
#'
#' @details
#' The function performs a join between sample and calibration tables via
#' `groupIndices`. Samples outside the linear range are classified as:
#' \itemize{
#'   \item BLOL: below lower limit of linearity
#'   \item ULOL: above upper limit of linearity
#' }
#' If calibration quality (`aboveR2`) is insufficient, the status is set to FALSE.
#'
#' @import data.table
#' @export
getLRstatus <- function(dats, datCal, y){

  data.table::setDT(dats)
  data.table::setDT(datCal)

  dat <- data.table::copy(dats)


  dat <- dat[datCal, on = "groupIndices"]
  if(any(startsWith(names(dat), "i."))){
  dat <- dat[, grep("^i\\.", names(dat), value=TRUE) := NULL]
  }
  dat$Status_LR <- NA
  dat$Status_LR = data.table::between(lower = dat$RangeStartY, x = dat[[y]], upper = dat$RangeEndY, NAbounds = NA)

  dat$Status_LR[dat[[y]] < dat$RangeStartY] <- "BLOL"
  dat$Status_LR[dat[[y]] > dat$RangeEndY] <- "ULOL"


  dat$Status_LR[dat$aboveR2 %in% c(NA,FALSE)] <- FALSE

  dats <- dats[dat[ ,.( IDintern,LRFlag, Status_LR, RangeStartY, RangeEndY)], on = "IDintern"]
  if(any(startsWith(names(dats), "i."))){
  dats <- dats[, grep("^i\\.", names(dats), value=TRUE) := NULL]
}
  return(dats)

}



#' Estimate concentrations from calibration curve parameters
#'
#' @description
#' Back-calculates sample concentrations using calibration curve parameters
#' (linear regression: y = intercept + slope * x).
#'
#' Depending on whether log-transformation is assumed (`INVERSE_Y`),
#' concentrations are either:
#' - back-transformed from log-space, or
#' - directly calculated in linear space.
#'
#' Finally, concentrations are scaled relative to a reference standard.
#'
#' @param dats Data table containing biological samples.
#' Must include:
#' \itemize{
#'   \item groupIndices (link to calibration curve)
#'   \item IDintern (sample identifier)
#' }
#'
#' @param datCal Data table with calibration parameters.
#' Must include:
#' \itemize{
#'   \item groupIndices
#'   \item Intercept (regression intercept)
#'   \item slope (regression slope)
#'   \item sigma (residual standard deviation, if log model is used)
#'   \item Compound (feature identifier)
#' }
#'
#' @param y Character. Name of response variable used in calibration model.
#'
#' @param INVERSE_Y Logical or NULL.
#' If not NULL, assumes log-normal model and applies exponential back-transformation.
#'
#' @param NAME_Standard Character vector defining standard samples used for scaling.
#'
#' @param COL_expConc Character. Column name containing expected concentrations.
#'
#' @return
#' Data table with added column:
#' \itemize{
#'   \item ConcentrationLR: back-calculated and scaled concentration estimate
#' }
#'
#' @details
#' Model assumptions:
#' \itemize{
#'   \item Linear model: y = intercept + slope * x
#'   \item Log model: x estimated via exponential back-transformation
#' }
#'
#' Final concentration is scaled using the minimum expected standard concentration:
#' \deqn{C_{scaled} = C \cdot \frac{1}{3} \cdot C_{min, standard}}
#'
#' @importFrom data.table setDT
#' @export
getConc <- function(dats, datCal, y, INVERSE_Y, NAME_Standard, COL_expConc){

  setDT(dats)
  setDT(datCal)

  dat <- data.table::copy(dats)


  dat <- dat[datCal, on = intersect(colnames(dat), colnames(datCal))]

  #y = b0 + b1*x

  if(!is.null(INVERSE_Y) ){

    dat[, ConcentrationLR := exp((get(y) - Intercept - sigma^2/2)/slope), by = .I]

    #dat$ConcentrationLR <- sapply(paste0(INVERSE_Y,"(",dat$ConcentrationLR,")"),function(i) eval(parse(text = i)))
  } else{

    dat[, ConcentrationLR := (get(y) - Intercept)/slope, by = .I]


  }
  dat[, ConcentrationLR := ConcentrationLR * 1/3 * min(as.numeric(get(COL_expConc))[Sample.Type %in% NAME_Standard], na.rm = TRUE), by = Compound]
  #dat$ConcentrationLR <- dat$ConcentrationLR/(3/dat$xfactor


  #dats <- dats[dat[ ,.(IDintern, ConcentrationLR)], on = "IDintern"]
  return(dat)

}


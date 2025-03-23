#' Title
#'
#' @param dats data table with information about the biological samples.
#' Required columns are IDintern, y, groupIndices.It is necessary
#' that the columns in dats and datCal have the same column names and that the
#' column groupIndices is compatible between samples and dilution/concentration curve data.
#' @param datCal data table with dilution/concentration curve information of the linearity.
#' Required columns are groupIndices, Intercept and slope. It is necessary
#' that the columns in dats and datCal have the same column names and that the
#' column groupIndices is compatible between samples and dilution/concentration curve data.
#' @param y
#'
#' @return
#' @export
#'
#' @examples
getLRstatus <- function(dats, datCal, y){

  data.table::setDT(dats)
  data.table::setDT(datCal)

  dat <- data.table::copy(dats)


  dat <- dat[datCal, on = "groupIndices"]
  if(any(startsWith(names(dat), "i."))){
  dat <- dat[, grep("^i\\.", names(dat), value=TRUE) := NULL]
  }
  dat$Status_LR <- NA
  dat$Status_LR = data.table::between(lower = dat$LRStartY, x = dat[[y]], upper = dat$LREndY, NAbounds = NA)

  dat$Status_LR[dat[[y]] < dat$LRStartY] <- "BLOL"
  dat$Status_LR[dat[[y]] > dat$LREndY] <- "ULOL"


  dat$Status_LR[dat$aboveR2 %in% c(NA,FALSE)] <- FALSE

  dats <- dats[dat[ ,.( IDintern,LRFlag, Status_LR, LRStartY, LREndY)], on = "IDintern"]
  if(any(startsWith(names(dats), "i."))){
  dats <- dats[, grep("^i\\.", names(dats), value=TRUE) := NULL]
}
  return(dats)

}




#' Title
#'
#' @param dats data table with information about the biological samples.
#' Required columns are IDintern, y, groupIndices.It is necessary
#' that the columns in dats and datCal have the same column names and that the
#' column groupIndices is compatible between samples and dilution/concentration curve data.
#' @param datCal data table with dilution/concentration curve information of the linearity.
#' Required columns are groupIndices, Intercept and slope. It is necessary
#' that the columns in dats and datCal have the same column names and that the
#' column groupIndices is compatible between samples and dilution/concentration curve data.
#' @param y
#' @param INVERSE_Y
#'
#' @return
#' @export
#'
#' @examples
getConc <- function(dats, datCal, y, INVERSE_Y, NAME_Standard, COL_expConc){

  setDT(dats)
  setDT(datCal)

  dat <- data.table::copy(dats)


  dat <- dat[datCal, on = intersect(colnames(dat), colnames(datCal))]

  #y = b0 + b1*x
  dat[, ConcentrationLR := (get(y) - Intercept)/slope, by = .I]

  if(!is.null(INVERSE_Y) ){
    dat$ConcentrationLR <- sapply(paste0(INVERSE_Y,"(",dat$ConcentrationLR,")"),function(i) eval(parse(text = i)))
  }
  dat[, ConcentrationLR := ConcentrationLR/(3/min(as.numeric(get(COL_expConc))[Sample.Type %in% NAME_Standard], na.rm = TRUE)), by = Compound]
  #dat$ConcentrationLR <- dat$ConcentrationLR/(3/dat$xfactor


  #dats <- dats[dat[ ,.(IDintern, ConcentrationLR)], on = "IDintern"]
  return(dat)

}


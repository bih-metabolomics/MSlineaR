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

  setDT(dats)
  setDT(datCal)

  dat <- data.table::copy(dats)

  dat <- dat[datCal, on = "groupIndices"]
  if(any(startsWith(names(dat), "i."))){
  dat <- dat[, grep("^i\\.", names(dat), value=TRUE) := NULL]
}
  dat$Status_LR = data.table::between(lower = dat$LRStartY, x = dat[[y]], upper = dat$LREndY, NAbounds = NA)


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
getConc <- function(dats, datCal, y, INVERSE_Y, x){

  setDT(dats)
  setDT(datCal)

  dat <- data.table::copy(dats)


  dat <- dat[datCal, on = "groupIndices"]

  #y = b0 + b1*x
  dat[, ConcentrationLR := (get(y) - Intercept)/slope, by = .I]

  if(!is.na(INVERSE_Y) & INVERSE_Y !=""){
    dat$ConcentrationLR <- sapply(paste0(INVERSE_Y,"(",dat$ConcentrationLR,")"),function(i) eval(parse(text = i)))
  }
  dat$ConcentrationLR <- dat$ConcentrationLR - get(x)

  dats <- dats[dat[ ,.(IDintern, ConcentrationLR)], on = "IDintern"]
  return(dats)

}


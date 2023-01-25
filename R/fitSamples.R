#' Title
#'
#' @param dats samples compare to the dilution series
#' @param datCal output object from AssessLinearity function$summaryFDS
#' @param COLNAMES
#'
#' @return
#' @export
#'
#' @examples
getLRstatus <- function(dats, datCal, COLNAMES = c(ID = "Compound", Replicate = "Batch", Y = "Area")){

  setDT(dats)
  setDT(datCal)

  dats$uniqueID <- 1:nrow(dats)
  dat <- data.table::copy(dats)
  data.table::setnames(dat,old =  COLNAMES[1:2], new = c(names(datCal)[2], names(datCal)[3]))

  dat <- dat[datCal, on = names(datCal)[2:3]]

  dat$Status_LR = data.table::between(lower = dat$LRStartY, x = dat[, get(COLNAMES[["Y"]])], upper = dat$LREndY)

  data.table::setnames(dat,new =  COLNAMES[1:2], old = c(names(datCal)[2], names(datCal)[3]))

  dats <- dats[dat[ ,.(uniqueID, LRFlag, Status_LR)], on = "uniqueID"]
  return(dats[, uniqueID := NULL])

}




#' Title
#'
#' @param dats samples compare to the dilution series
#' @param datCal [output object from AssessLinearity function]$summaryFDS
#' @param COLNAMES
#'
#' @return
#' @export
#'
#' @examples
getConc <- function(dats, datCal, COLNAMES = c(ID = "Compound", Replicate = "Batch", Y = "Area")){

  setDT(dats)
  setDT(datCal)

  dats$uniqueID <- 1:nrow(dats)
  dat <- data.table::copy(dats)
  data.table::setnames(dat,old =  COLNAMES[1:2], new = c(names(datCal)[2], names(datCal)[3]))

  dat <- dat[datCal, on = names(datCal)[2:3]]

  #y = b0 + b1*x
  dat[, ConcentrationLR := exp((log(get(COLNAMES[["Y"]])) - Intercept)/slope), by = .I]

  data.table::setnames(dat,new =  COLNAMES[1:2], old = c(names(datCal)[2], names(datCal)[3]))

  dats <- dats[dat[ ,.(uniqueID, ConcentrationLR)], on = "uniqueID"]
  return(dats[, uniqueID := NULL])

}


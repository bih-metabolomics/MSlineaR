#' Prepare raw data
#'
#' @param dat Data frame with raw data. Mandatory columns are:
#' - ID,
#' - IDintern
#' - Batch,
#' - Replicate,
#' - MZ
#' - RT
#' - Concentration
#' - Intensity
#'
#'
#' @return Data frame with columns:
#' - IDintern,
#' - ID,
#' - IntensityNorm,
#' - DP,
#' - Comment,
#' - pch,
#' - color,
#' - ConcentrationRaw,
#' - IntensityRaw,
#' - ConcentrationLog,
#' - IntensityLog,
#' - IntensityRaw
#' - DilutionPoint
#'
#' @export
#'
#' @examples
#'
checkData <- function(dat, MIN_FEATURE = 3, LOG_TRANSFORM = TRUE, nCORE = 1, ...){

  if (!"REPLICATE" %in% names(COLNAMES)) {
    COLNAMES["REPLICATE"] <- "REPLICATE"
    dat[, REPLICATE := "1"]
  }

  data.table::set(dat, j = COLNAMES[["ID"]], value = as.character(dat[[COLNAMES[["ID"]]]]))
  data.table::set(dat, j = COLNAMES[["REPLICATE"]], value = as.character(dat[[COLNAMES[["REPLICATE"]]]]))
  data.table::set(dat, j = COLNAMES[["X"]], value = as.numeric(dat[[COLNAMES[["X"]]]]))
  data.table::set(dat, j = COLNAMES[["Y"]], value = as.numeric(dat[[COLNAMES[["Y"]]]]))


  data.table::setorderv(dat, c(COLNAMES[["ID"]], COLNAMES[["REPLICATE"]], COLNAMES[["X"]], COLNAMES[["Y"]]))
  dat[ , IDintern := paste0("s", 1:.N)]
  dat <- dat[ , c( "IDintern", COLNAMES[["ID"]], COLNAMES[["REPLICATE"]], COLNAMES[["X"]], COLNAMES[["Y"]]), with = F]

  # tests for Input parameter

  stopifnot(exprs = {
    "missing columns" = all(COLNAMES %in% colnames(dat))
    "all values of X and Y need to be from type double and positive." = is.double(dat[[COLNAMES[["X"]]]]) & all(dat[[COLNAMES[["X"]]]] >= 0)
    "all values of X and Y need to be from type double and positive." = is.double(dat[[COLNAMES[["Y"]]]]) & all(dat[[COLNAMES[["Y"]]]] >= 0)
    "Column ID needs to be from type character " = is.character(dat[[COLNAMES[["ID"]]]])
    "Column REPLICATE needs to be from type character " = is.character(dat[[COLNAMES[["REPLICATE"]]]])
    "Parameter LOG_TRANSFORM need to be from type logical" = is.logical(LOG_TRANSFORM)
    "Parameter nCore needs to be positive" = nCORE > 0
    "Parameter MIN_FEATURE needs to be greater or equal than 3" = MIN_FEATURE >= 3

  })



  return(dat)
}


#' Title
#'normalizing, centralizing, log transforming
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
prepareData <- function(dat){

  stopifnot(exprs = {
    "data needs 5 columns called: IDintern, ID, REPLICATE, X, Y" = dim(dat)[2] == 5
    })

  processed <- data.table::copy(dat)
  # rename
  data.table::setnames(x = processed, old = colnames(processed), new = c( "IDintern", "ID", "REPLICATE", "X", "Y"))

  setorderv(processed, c("ID", "REPLICATE", "X"))

  processed[, ":="(Comment = NA, pch = fcase(!is.na(Y), 19), color = fcase(!is.na(Y), "black"))]
  processed[, ":="(YNorm = Y / max(Y, na.rm = T) * 100,
             YLog = log(Y),
             XLog = log(X),
             DilutionPoint = 1:.N,
             groupIndices = .GRP), by = c("ID", "REPLICATE")]
  processed[ , c( "IDintern","ID", "REPLICATE", "X", "Y", "Comment", "YNorm", "YLog", "XLog", "DilutionPoint", "groupIndices", "pch", "color"), with = F]

  return(processed)
}


# progressbar
pbapply::pboptions(type = "timer", char = "[=-]", style = 5)

progressr::handlers(global = TRUE)
progressr::handlers(
  progressr::handler_progress(
    format   = ":current/:total [:bar] :percent in :elapsed ETA: :eta",
    width    = 60,
    complete = "+"
  )
)

# if(nCORE > 1){


# handlers(global = TRUE)
my_fcn <- function(xs, func, inputData, ...) {
  p <- progressr::progressor(along = xs)
  y <- furrr::future_map(xs, function(x) {
    p(sprintf("x=%g", x))
    func(setDT(inputData)[groupIndices %in% unique(groupIndices)[x]], ...)
    # groupInd <- unique(inputData$groupIndices)[x]
    # inputData = filter(inputData, groupIndices %in% groupInd)
    # inputData[groupIndices %in% x, func(.SD, ...)]
    # func(inputData, ...)
  }, .progress = F)
}

# } else{
my_fcn2 <- function(xs, func, inputData, ...) {
  p <- progressr::progressor(along = xs)
  y <- map(xs, function(x) {
    p(sprintf("x=%g", x))
    # inputData[groupIndices %in% x, func(.SD, ...)]
    func(setDT(inputData)[groupIndices %in% unique(groupIndices)[x]], ...)
    # groupInd <- unique(inputData$groupIndices)[x]
    # inputData = filter(inputData, groupIndices %in% groupInd)
    # func(inputData, ...)
  })
}
# }

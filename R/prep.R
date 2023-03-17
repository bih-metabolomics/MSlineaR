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
checkData <- function(DAT, ...){

  dat <- DAT
  data.table::setDT(dat)

  stopifnot(exprs = {
    "'input_dat' has less rows than 'min_feature' required" = dim(dat)[1] >= MIN_FEATURE
    "'input_dat' needs minimum 4 columns defining the ID, the x and y values and the sample type of the data" = dim(dat)[2] >= 4

  })

  if (!"Batch" %in% names(COLNAMES)) {
    COLNAMES["Batch"] <- "Batch"
    dat[, Batch := "1"]
  }

  data.table::set(dat, j = COLNAMES[["ID"]], value = as.character(dat[[COLNAMES[["ID"]]]]))
  data.table::set(dat, j = COLNAMES[["Batch"]], value = as.character(dat[[COLNAMES[["Batch"]]]]))
  data.table::set(dat, j = COLNAMES[["X"]], value = as.numeric(dat[[COLNAMES[["X"]]]]))
  data.table::set(dat, j = COLNAMES[["Y"]], value = as.numeric(dat[[COLNAMES[["Y"]]]]))
  data.table::set(dat, j = COLNAMES[["Sample_type"]], value = as.character(dat[[COLNAMES[["Sample_type"]]]]))


  data.table::setorderv(dat, c(COLNAMES[["ID"]], COLNAMES[["Sample_type"]], COLNAMES[["Batch"]], COLNAMES[["X"]], COLNAMES[["Y"]]))
  dat[ , IDintern := paste0("s", 1:.N)]
  dat <- dat[ , c( "IDintern", COLNAMES[["ID"]], COLNAMES[["Sample_type"]], COLNAMES[["Batch"]], COLNAMES[["X"]], COLNAMES[["Y"]]), with = F]

  # tests for Input parameter

  stopifnot(exprs = {

    "missing columns" = all(COLNAMES %in% colnames(dat))
    "all values of 'column_X' and 'column_Y' need to be from type double and positive." =
      is.double(dat[[COLNAMES[["X"]]]]) & all(dat[[COLNAMES[["X"]]]] >= 0, na.rm = T)
    "all values of 'column_X' and 'column_Y' need to be from type double and positive." =
      is.double(dat[[COLNAMES[["Y"]]]]) & all(dat[[COLNAMES[["Y"]]]] >= 0, na.rm = T)
    "'column_ID' needs to be from type character " = is.character(dat[[COLNAMES[["ID"]]]])
    "'column_Batch' needs to be from type character " = is.character(dat[[COLNAMES[["Batch"]]]])
    "Argument 'transform' needs to be from type logical" = is.logical(TRANSFORM)
    "Argument 'first_outlier_detection' needs to be from type logical" = is.logical(FOD)
    "Argument 'second_outlier_detection' needs to be from type logical" = is.logical(SOD)
    "Argument 'trimming' needs to be from type logical" = is.logical(TRIMM)
    "Argument 'get_output' needs to be from type logical" = is.logical(GET_OUTPUT)
    "Argument 'calculate_concentration' needs to be from type logical" = is.logical(CAL_CONC)
    "Argument 'get_linearity_status_samples' needs to be from type logical" = is.logical(GET_LR_STATUS)
    "Argument 'nCore' needs to be positive" = nCORE > 0
    "Argument 'min_feature' needs to be greater or equal than 3" = MIN_FEATURE >= 3
    "Argument 'LR_sd_res_factor' needs to be from type integer and positive" = is.wholenumber(LR_SD_RES_FACTOR) & LR_SD_RES_FACTOR >= 0
    "Argument 'sample_type_serial' was not found in column 'column_sample_type'" =
        any(dat[[column_sample_type]] %in% CALIBRANTS)
    "Argument 'R2_min' needs to be from type double and in the range between 0 and 1" =
      data.table::between(x = R2_MIN, lower = 0, upper = 1)

  })

  if(TRANSFORM %in% TRUE){
    stopifnot(exprs = {
      "Argument 'transform_X' and 'transform_Y' need to be a String or NA, e.g. 'log10'" =
        (is.character(TRANSFORM_X)| is.na(TRANSFORM_X)) & (is.character(TRANSFORM_Y)| is.na(TRANSFORM_Y))

    })
  }

  if(FOD %in% TRUE){
    stopifnot(exprs = {
      "For the Argument FOD_model only one or a combination of c('linear', 'logistic', 'quadratic') are allowed. " =
        FOD_MODEL %in% c('linear', 'logistic', 'quadratic')
      "Argument 'FOD_sdres_min' and 'FOD_stdres_max' need to be from type Integer and positive" =
        is.wholenumber(FOD_SDRES_MIN) & FOD_SDRES_MIN >= 0  & is.wholenumber(FOD_STDRES_MAX) & FOD_STDRES_MAX >= 0

    })
  }


  if(SOD %in% TRUE){
    stopifnot(exprs = {
      "For the Argument SOD_model only one or a combination of c('linear', 'logistic', 'quadratic') are allowed. " =
        SOD_MODEL %in% c('linear', 'logistic', 'quadratic')
      "Argument 'SOD_sdres_min' and 'SOD_stdres_max' need to be from type Integer and positive" =
        is.wholenumber(SOD_SDRES_MIN) & SOD_SDRES_MIN >= 0  & is.wholenumber(SOD_STDRES_MAX) & SOD_STDRES_MAX >= 0

    })
  }

if(get_linearity_status_samples %in% TRUE){

  stopifnot(exprs = {
    "Argument 'sample_type_sample' was not found in column 'column_sample_type'" =
      any(dat[[column_sample_type]] %in% SAMPLE)

  })
}




  return(dat)
}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#' Title
#'normalizing, centralizing, log transforming
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
prepareData <- function(dat,...){

  stopifnot(exprs = {
    "data need to have 6 columns" = dim(dat)[2] == 6
    })
data.table::setDT(dat)
  processed <- data.table::copy(dat)
  # rename
  data.table::setnames(x = processed, old = colnames(processed), new = c( "IDintern", "ID","Sample.Type", "Batch", "X", "Y"))

  data.table::setorderv(processed, c("ID", "Batch", "X"))

  processed[, ":="(Comment = NA, pch = data.table::fcase(!is.na(Y), 19), color = data.table::fcase(!is.na(Y), "black"))]
  processed[, ":="(#YNorm = Y / max(Y, na.rm = T) * 100,
             #Y_trans = log(Y),
             #X_trans = log(X),
             DilutionPoint = 1:.N,
             groupIndices = .GRP), by = c("ID", "Batch")]

  if(TRANSFORM %in% TRUE & !is.na(TRANSFORM_X)){
    processed$X_trans = get(TRANSFORM_X)(processed$X)
    processed$X_trans[is.infinite(processed$X_trans)] <- NA
  }

  if(TRANSFORM %in% TRUE & !is.na(TRANSFORM_Y)){
    processed$Y_trans = get(TRANSFORM_Y)(processed$Y)
    processed$Y_trans[is.infinite(processed$Y_trans)] <- NA
  }

  #processed <- processed[ , c("groupIndices", "IDintern","ID", "Batch", "X", "Y", "Comment", "YLog", "XLog", "DilutionPoint", "pch", "color"), with = F]

  return(processed)
}



# if(nCORE > 1){


# handlers(global = TRUE)
#' Title
#'
#' @param xs
#' @param func
#' @param inputData
#' @param ...
#'
#' @return
#' @import foreach
#' @export
#'
#' @examples
#'

my_fcn <- function(xs, func, inputData, ...) {
  #parallel::clusterExport(cl, exportObjects)
  cl <- parallel::makePSOCKcluster(nCORE)
   doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  #on.exit(closeAllConnections())
  #on.exit(parallel::stopCluster(cl))
  pb <- txtProgressBar(min=1, max=max(xs), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
   #pb <- progressr::progressor(along = xs)
   #progress <- function(i)  pb(sprintf("x=%g", i))
  opts <- list(progress=progress)
  y <- foreach::foreach(i = xs,   .options.snow=opts) %dopar%
    {func(data.table::setDT(inputData)[inputData$groupIndices %in% unique(inputData$groupIndices)[i]], ...)}
  close(pb)
  #parallel::stopCluster(cl)
  #parallel::stopCluster(cl)
  return(y)
  }

my_fcn6 <- function(xs, func, inputData, ...) {
  #parallel::clusterExport(cl, exportObjects)
  #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
  #on.exit(parallel::stopCluster(cl))
  #doSNOW::registerDoSNOW(cl)
#  pb <- progressr::progressor(along = xs)
#  progress <- function(i)  pb(sprintf("x=%g", i))
 # opts <- list(progress=progress)
  y <- foreach::foreach(i = xs) %do% {func(data.table::setDT(inputData)[inputData$groupIndices %in% unique(inputData$groupIndices)[i]], ...)}
  #parallel::stopCluster(cl)
  #return(y)
}




# my_fcn <- function(xs, func, inputData, ...) {
#   p <- progressr::progressor(along = xs)
#   y <- furrr::future_map(xs, function(i) {
#     p(sprintf("x=%g", i))
#     func(data.table::setDT(inputData)[inputData$groupIndices %in% unique(inputData$groupIndices)[i]], ...)
#     # groupInd <- unique(inputData$groupIndices)[x]
#     # inputData = filter(inputData, groupIndices %in% groupInd)
#     # inputData[groupIndices %in% x, func(.SD, ...)]
#     # func(inputData, ...)
#   }, .progress = F)
# }

# } else{
my_fcn2 <- function(xs, func, inputData, ...) {
  p <- progressr::progressor(along = xs)
  y <- map(xs, function(i) {
    p(sprintf("x=%g", i))
    # inputData[groupIndices %in% x, func(.SD, ...)]
    func(data.table::setDT(inputData)[groupIndices %in% unique(groupIndices)[i]], ...)
    # groupInd <- unique(inputData$groupIndices)[x]
    # inputData = filter(inputData, groupIndices %in% groupInd)
    # func(inputData, ...)
  })
}
# }

#' Prepare raw data
#'
#' @param MIN_FEATURE
#' @param TYPE
#' @param QC
#' @param QC_REF
#' @param BLANK
#' @param SAMPLE
#' @param SAMPLE_ID
#' @param CALIBRANTS
#' @param COLNAMES
#' @param Y_SAMPLE
#' @param TRANSFORM
#' @param TRANSFORM_X
#' @param INVERSE_X
#' @param TRANSFORM_Y
#' @param INVERSE_Y
#' @param FOD
#' @param FOD_MODEL
#' @param FOD_SDRES_MIN
#' @param FOD_STDRES_MAX
#' @param TRIMM
#' @param SOD
#' @param SOD_MODEL
#' @param SOD_SDRES_MIN
#' @param SOD_STDRES_MAX
#' @param LR_SD_RES_FACTOR
#' @param R2_MIN
#' @param BATCH_HARMONIZATION
#' @param CAL_CONC
#' @param GET_LR_STATUS
#' @param nCORE
#' @param GET_OUTPUT
#' @param PREFIX
#' @param OUTPUT_DIR
#' @param dat
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
checkData <- function(dat, MIN_FEATURE = parent.frame()$MIN_FEATURE, TYPE = parent.frame()$TYPE, QC = parent.frame()$QC,
                      BLANK = parent.frame()$BLANK,
                      SAMPLE = parent.frame()$SAMPLE,Y_SAMPLE = parent.frame()$Y_SAMPLE,
                      CALIBRANTS = parent.frame()$CALIBRANTS, COLNAMES = parent.frame()$COLNAMES,
                      DILUTION_FACTOR = parent.frame()$DILUTION_FACTOR,
                      NOISE = parent.frame()$NOISE,
                      TRANSFORM = parent.frame()$TRANSFORM,TRANSFORM_X = parent.frame()$TRANSFORM_X, INVERSE_X = parent.frame()$INVERSE_X,
                      TRANSFORM_Y = parent.frame()$TRANSFORM_Y,INVERSE_Y = parent.frame()$INVERSE_Y,
                      FOD = parent.frame()$FOD,FOD_MODEL = parent.frame()$FOD_MODEL, FOD_SDRES_MIN = parent.frame()$FOD_SDRES_MIN, FOD_STDRES_MAX = parent.frame()$FOD_STDRES_MAX,
                      TRIMM = parent.frame()$TRIMM,
                      SOD = parent.frame()$SOD,SOD_MODEL = parent.frame()$SOD_MODEL, SOD_SDRES_MIN = parent.frame()$SOD_SDRES_MIN, SOD_STDRES_MAX = parent.frame()$SOD_STDRES_MAX,
                      LR_SD_RES_FACTOR = parent.frame()$LR_SD_RES_FACTOR, R2_MIN = parent.frame()$R2_MIN,
                      BATCH_HARMONIZATION = parent.frame()$BATCH_HARMONIZATION,
                      CAL_CONC = parent.frame()$CAL_CONC,GET_LR_STATUS = parent.frame()$GET_LR_STATUS,
                      nCORE = parent.frame()$nCORE,
                      GET_OUTPUT = parent.frame()$GET_OUTPUT, IMG_OUTPUT_DIR = parent.frame()$IMG_OUTPUT_DIR){

  data.table::setDT(dat)


  if(dim(dat)[1] < MIN_FEATURE) rlang::abort("'input_dat' has less rows than argument 'min_feature' required")
  if(dim(dat)[2] < 4) rlang::abort("'input_dat' needs minimum 4 columns defining the ID, the x and y values and the sample type of the data")

  if (!"Batch" %in% names(COLNAMES)) {
    COLNAMES["Batch"] <- "Batch"
    dat[, Batch := "1"]
  } else if("Batch" %in% names(COLNAMES) & is.na(COLNAMES[["Batch"]])){
    COLNAMES["Batch"] <- "Batch"
    dat[, Batch := "1"]
  }


  if (!"Class" %in% names(COLNAMES)) {
    COLNAMES["Class"] <- "Class"
    dat[, Class := "1"]
  } else if("Class" %in% names(COLNAMES) & is.na(COLNAMES[["Class"]])){
    COLNAMES["Class"] <- "Class"
    dat[, Class := "1"]
  }

  if(!all(COLNAMES %in% colnames(dat)))rlang::abort("missing columns")


  data.table::set(dat, j = COLNAMES[["ID"]], value = as.character(dat[[COLNAMES[["ID"]]]]))
  data.table::set(dat, j = COLNAMES[["Batch"]], value = as.character(dat[[COLNAMES[["Batch"]]]]))
  data.table::set(dat, j = COLNAMES[["Class"]], value = as.character(dat[[COLNAMES[["Class"]]]]))
  data.table::set(dat, j = COLNAMES[["X"]], value = as.numeric(dat[[COLNAMES[["X"]]]]))
  data.table::set(dat, j = COLNAMES[["Y"]], value = as.numeric(dat[[COLNAMES[["Y"]]]]))
  data.table::set(dat, j = COLNAMES[["Sample_type"]], value = as.character(dat[[COLNAMES[["Sample_type"]]]]))


  data.table::setorderv(dat, c(COLNAMES[["ID"]], COLNAMES[["Sample_ID"]], COLNAMES[["Sample_type"]],COLNAMES[["Class"]], COLNAMES[["Batch"]], COLNAMES[["X"]], COLNAMES[["Y"]]))
  dat[ , IDintern := paste0("s", 1:.N)]
  dat[ , groupIndices := .GRP ,by = c(COLNAMES[["ID"]], COLNAMES[["Batch"]])]
  dat <- dat[ , c( "IDintern", "groupIndices", COLNAMES[["ID"]],  COLNAMES[["Sample_ID"]], COLNAMES[["Sample_type"]],COLNAMES[["Class"]], COLNAMES[["Batch"]], COLNAMES[["X"]], COLNAMES[["Y"]]), with = F]

  # tests for Input parameter

  if(!(is.double(dat[[COLNAMES[["X"]]]]) & all(dat[[COLNAMES[["X"]]]] >= 0, na.rm = T))) rlang::abort("all values of 'column_X' and 'column_Y' need to be from type double and positive.")
  if(!(is.double(dat[[COLNAMES[["X"]]]]) & all(dat[[COLNAMES[["Y"]]]] >= 0, na.rm = T))) rlang::abort("all values of 'column_X' and 'column_Y' need to be from type double and positive.")
  if(!is.character(dat[[COLNAMES[["ID"]]]])) rlang::abort("'column_ID' needs to be from type character ")
  if(!is.character(dat[[COLNAMES[["Batch"]]]])) rlang::abort("'column_Batch' needs to be from type character ")
  if(!is.character(dat[[COLNAMES[["Class"]]]])) rlang::abort("'column_Class' needs to be from type character ")

  if(!is.logical(FOD)| is.na(FOD)) rlang::abort("Argument 'first_outlier_detection' needs to be from type logical")
  if(FOD %in% TRUE) if(!all(FOD_MODEL %in% c('linear', 'logistic', 'quadratic')))rlang::abort("For the Argument 'FOD_model' only one or a combination of 'linear', 'logistic' or 'quadratic' are allowed.")
  if(FOD %in% TRUE & (!is.numeric(FOD_SDRES_MIN) | !is.numeric(FOD_STDRES_MAX)))  rlang::abort("Argument 'FOD_sdres_min' and 'FOD_stdres_max' need to be from type Integer and positive")
  if(FOD %in% TRUE) if(!(is.wholenumber(FOD_SDRES_MIN) & FOD_SDRES_MIN >= 0  & is.wholenumber(FOD_STDRES_MAX) & FOD_STDRES_MAX >= 0)) rlang::abort("Argument 'FOD_sdres_min' and 'FOD_stdres_max' need to be from type Integer and positive")

  if(!is.logical(TRIMM) | is.na(TRIMM)) rlang::abort("Argument 'trimming' needs to be from type logical")

  if(!is.logical(SOD) | is.na(SOD)) rlang::abort("Argument 'second_outlier_detection' needs to be from type logical")
  if(SOD %in% TRUE) if(! all(SOD_MODEL %in% c('linear', 'logistic', 'quadratic')))rlang::abort("For the Argument SOD_model only one or a combination of 'linear', 'logistic', 'quadratic' are allowed. ")
  if(SOD %in% TRUE & (!is.numeric(SOD_SDRES_MIN) | !is.numeric(SOD_STDRES_MAX))) rlang::abort ("Argument 'SOD_sdres_min' and 'SOD_stdres_max' need to be from type Integer and positive")
  if(SOD %in% TRUE) if(!(is.wholenumber(SOD_SDRES_MIN) & SOD_SDRES_MIN >= 0  & is.wholenumber(SOD_STDRES_MAX) & SOD_STDRES_MAX >= 0)) rlang::abort ("Argument 'SOD_sdres_min' and 'SOD_stdres_max' need to be from type Integer and positive")

  if(!is.logical(GET_OUTPUT) | is.na(GET_OUTPUT)) rlang::abort("Argument 'get_output' needs to be from type logical")
  if(!is.logical(CAL_CONC) | is.na(CAL_CONC)) rlang::abort("Argument 'calculate_concentration' needs to be from type logical")
  if(!is.logical(GET_LR_STATUS) | is.na(GET_LR_STATUS)) rlang::abort("Argument 'get_linearity_status_samples' needs to be from type logical")
  if(!nCORE > 0) rlang::abort("Argument 'nCore' needs to be positive")
  if(!MIN_FEATURE >= 3) rlang::abort("Argument 'min_feature' needs to be greater or equal than 3")
  if(!is.numeric(LR_SD_RES_FACTOR)) rlang::abort("Argument 'LR_sd_res_factor' needs to be from type integer and positive")
  if(!is.wholenumber(LR_SD_RES_FACTOR) | LR_SD_RES_FACTOR < 0) rlang::abort("Argument 'LR_sd_res_factor' needs to be from type integer and positive")

  if(!any(dat[[COLNAMES[["Sample_type"]]]] %in% CALIBRANTS)) rlang::abort("Argument 'sample_type_serial' was not found in column 'column_sample_type'")
  if(!data.table::between(x = R2_MIN, lower = 0, upper = 1)) rlang::abort("Argument 'R2_min' needs to be from type double and in the range between 0 and 1")
  if(!TYPE %in% c("untargeted", "targeted")) rlang::abort("Argument 'analysis_type' need to be either 'untargeted' or 'targeted'")
  if(TYPE %in% "untargeted")if(!(!is.null(DILUTION_FACTOR) | DILUTION_FACTOR != " " | !is.na(DILUTION_FACTOR))) rlang::abort("Argument 'dilution_factor ' must be provided for untargeted analysis")

  if(!is.logical(TRANSFORM) | is.na(TRANSFORM)) rlang::abort("Argument 'transform' needs to be from type logical")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_X)) if(!is.character(TRANSFORM_X)) rlang::abort("Argument 'transform_X' and 'transform_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_X)) if(!class(get(TRANSFORM_X)) == "function") rlang::abort("Argument 'transform_X' and 'transform_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)) if(!is.character(TRANSFORM_Y)) rlang::abort("Argument 'transform_X' and 'transform_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)) if(!class(get(TRANSFORM_Y)) == "function") rlang::abort("Argument 'transform_X' and 'transform_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(INVERSE_X))   if(!is.character(INVERSE_X)) rlang::abort("Argument 'inverse_X' and 'inverse_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(INVERSE_X))   if(!class(get(INVERSE_X)) == "function") rlang::abort("Argument 'inverse_X' and 'inverse_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(INVERSE_Y))   if(!is.character(INVERSE_Y)) rlang::abort("Argument 'inverse_X' and 'inverse_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(INVERSE_Y))   if(!class(get(INVERSE_Y)) == "function") rlang::abort("Argument 'inverse_X' and 'inverse_Y' needs to be a String indicating a function, e.g. 'log10', disable with NULL")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_X)) if(is.null(INVERSE_X)) rlang::abort("Argument 'inverse_X' must be provided if Argument 'transform_X' is provided")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)) if(is.null(INVERSE_Y)) rlang::abort("Argument 'inverse_Y' must be provided if Argument 'transform_Y' is provided")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_X)) if(!as.integer(get(INVERSE_X)(get(TRANSFORM_X)(10))) == 10) rlang::abort( "Argument 'inverse_X' must be the reverse function of Argument 'transform_X'")
  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)) if(!as.integer(get(INVERSE_Y)(get(TRANSFORM_Y)(10))) == 10) rlang::abort( "Argument 'inverse_Y' must be the reverse function of Argument 'transform_Y'")


  if(!is.null(SAMPLE)) if(!any(dat[[COLNAMES[["Sample_type"]]]] %in% SAMPLE))rlang::abort("Argument 'sample_type_sample' was not found in column 'column_sample_type'")
  if(!is.null(SAMPLE)) if(!any(colnames(dat) %in% Y_SAMPLE))rlang::abort("Column 'column_Y_sample' was not found in input data.")

  if(!is.null(QC)) if(!all( QC %in% dat[[COLNAMES[["Sample_type"]]]])) rlang::abort("Argument 'sample_type_QC' was not found in column 'column_sample_type'")
  if(!is.numeric(DILUTION_FACTOR)) rlang::abort("Argument 'dilution_factor' needs to be from type numeric.")
  if(!is.logical(BATCH_HARMONIZATION) | is.na(BATCH_HARMONIZATION)) rlang::abort("Argument 'Batch_harmonization' needs to be from type logical")
  if(!is.null(NOISE)) if(!is.numeric(NOISE) | NOISE < 1) rlang::abort("Argument 'signal_blank_ration' needs to be from type integer and greater than 0.")
  if(!is.null(NOISE)) if(!BLANK %in% dat[[COLNAMES[["Sample_type"]]]]) rlang::abort("Argument 'sample_type_blank' was not found in column 'column_sample_type'")




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
prepareData <- function(dat,
                        TRANSFORM = parent.frame()$TRANSFORM,
                        TRANSFORM_X = parent.frame()$TRANSFORM_X,
                        TRANSFORM_Y = parent.frame()$TRANSFORM_Y,
                        DILUTION_FACTOR = parent.frame()$DILUTION_FACTOR,
                        COLNAMES = parent.frame()$COLNAMES,
                        TYPE = parent.frame()$TYPE){

  if(dim(dat)[2] != 9) rlang::abort("data need to have 9 columns, please use funtion 'checkData' before to check all necessary input arguments")



  data.table::setDT(dat)
  processed <- data.table::copy(dat)
  # rename
  data.table::setnames(x = processed, old = colnames(processed), new = c( "IDintern","groupIndices", "ID","Sample_ID", "Sample.Type","Class", "Batch",COLNAMES[["X"]], "Y"))
  data.table::setorderv(processed, c("ID", "Batch", COLNAMES[["X"]]))

  processed[ , ":="(Comment = NA, pch = data.table::fcase(!is.na(Y), 19), color = data.table::fcase(!is.na(Y), "black"))]
  processed[ , ":="(#YNorm = Y / max(Y, na.rm = T) * 100,
    #Y_trans = log(Y),
    #X_trans = log(X),
    DilutionPoint = 0 : (.N-1)),
    #groupIndices = .GRP),
    by = c("ID", "Batch")]

  processed$X = DILUTION_FACTOR^processed$DilutionPoint
  processed$DilutionPoint <- processed$DilutionPoint + 1

  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_X)){
    processed$X_trans = get(TRANSFORM_X)(processed$X)
    processed$X_trans[is.infinite(processed$X_trans)] <- NA
  }

  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)){
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

my_fcn <- function(nCORE, xs, func, inputData, ...) {
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

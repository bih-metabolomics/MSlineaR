#' Validate and standardise raw metabolomics input data
#'
#' @description
#' This function performs extensive validation and standardisation of raw
#' mass spectrometry feature tables prior to downstream analysis in MSlineaR.
#'
#' It ensures that required columns are present, data types are correct,
#' sample annotations are consistent, and all key parameters meet the
#' requirements of the linearity assessment workflow.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item Minimum number of features and required columns
#'   \item Presence and consistency of sample annotations
#'   \item Correct data types for numeric and categorical variables
#'   \item Validity of model and preprocessing parameters
#'   \item Consistency of transformation functions (if used)
#' }
#'
#' In addition, the function:
#' \itemize{
#'   \item Converts key identifiers to appropriate data types
#'   \item Creates a unique internal identifier (`IDintern`)
#'   \item Defines feature grouping variables (`groupIndices`)
#'   \item Ensures correct ordering of samples for downstream modelling
#' }
#'
#' This function should always be called before \code{prepareData()} or any
#' downstream MSlineaR analysis steps.
#'
#' @param dat A data.table or data.frame containing raw feature intensity data.
#'
#' @param MIN_FEATURE Integer. Minimum number of features required in the dataset.
#'
#' @param TYPE Character. Type of analysis: either \code{"untargeted"} or \code{"targeted"}.
#'
#' @param QC Character vector defining QC sample types.
#'
#' @param QC_REF Character vector defining reference QC samples (if applicable).
#'
#' @param BLANK Character vector defining blank sample types.
#'
#' @param SAMPLE Character vector defining biological sample types.
#'
#' @param SAMPLE_ID Character. Column name containing sample identifiers.
#'
#' @param CALIBRANTS Character vector defining calibration or dilution series samples.
#'
#' @param COLNAMES Named list defining required column mappings (e.g. Sample_ID,
#' Feature_ID, X, Y, Batch, Injection_order, Sample_type, Class).
#'
#' @param TRANSFORM Logical. Whether transformation of X/Y values should be applied.
#'
#' @param TRANSFORM_X Character. Name of transformation function applied to X (e.g. "log10").
#'
#' @param INVERSE_X Character. Inverse transformation function for X.
#'
#' @param TRANSFORM_Y Character. Name of transformation function applied to Y.
#'
#' @param INVERSE_Y Character. Inverse transformation function for Y.
#'
#' @param FOD Logical. Enable first outlier detection.
#'
#' @param FOD_MODEL Character vector specifying models used for first outlier detection
#' (allowed: "linear", "logistic", "quadratic").
#'
#' @param FOD_SDRES_MIN Numeric. Minimum standardised residual threshold for FOD.
#'
#' @param FOD_STDRES_MAX Numeric. Maximum standardised residual threshold for FOD.
#'
#' @param TRIMM Logical. Enable trimming of non-linear regions.
#'
#' @param SOD Logical. Enable second outlier detection.
#'
#' @param SOD_MODEL Character vector specifying models used for second outlier detection.
#'
#' @param SOD_SDRES_MIN Numeric. Minimum residual threshold for SOD.
#'
#' @param SOD_STDRES_MAX Numeric. Maximum residual threshold for SOD.
#'
#' @param LR_SD_RES_FACTOR Numeric. Scaling factor for defining linear range boundaries.
#'
#' @param R2_MIN Numeric (0–1). Minimum R² required for acceptance of a linear range.
#'
#' @param BATCH_HARMONIZATION Logical. If TRUE, features must pass criteria in all batches.
#'
#' @param GET_LR_STATUS Logical. If TRUE, returns sample-level linearity classification.
#'
#' @param nCORE Integer. Number of CPU cores for parallel processing.
#'
#' @param GET_OUTPUT Logical. If TRUE, additional diagnostic outputs are returned.
#'
#' @param IMG_OUTPUT_DIR Character. Directory for saving diagnostic plots.
#'
#' @return A validated and standardised data.table with harmonised column structure,
#' including:
#' \itemize{
#'   \item IDintern – unique sample-feature-batch identifier
#'   \item groupIndices – feature-batch grouping index
#'   \item Feature_ID – feature identifier
#'   \item Sample_ID – sample identifier
#'   \item Sample_type – sample annotation
#'   \item Class – sample class label
#'   \item Batch – batch identifier
#'   \item Injection_order – acquisition order
#'   \item X – concentration / dilution variable
#'   \item Y – intensity variable
#' }
#'
#' @examples
#' \dontrun{
#' dat_checked <- checkData(dat, MIN_FEATURE = 100, TYPE = "untargeted",
#'                          QC = c("QC"), BLANK = c("blank"),
#'                          CALIBRANTS = c("cal1","cal2"),
#'                          COLNAMES = list(Sample_ID="Sample",
#'                                          Feature_ID="Feature",
#'                                          X="Concentration",
#'                                          Y="Intensity",
#'                                          Batch="Batch",
#'                                          Injection_order="Order",
#'                                          Sample_type="Type",
#'                                          Class="Class"))
#' }
#' @keywords internal
#'
#' @import data.table
#' @importFrom rlang abort
checkData <- function(dat, MIN_FEATURE = parent.frame()$MIN_FEATURE, TYPE = parent.frame()$TYPE, QC = parent.frame()$QC,
                      BLANK = parent.frame()$BLANK,
                      SAMPLE = parent.frame()$SAMPLE,#Y_SAMPLE = parent.frame()$Y_SAMPLE,
                      CALIBRANTS = parent.frame()$CALIBRANTS, COLNAMES = parent.frame()$COLNAMES,

                      NOISE = parent.frame()$NOISE,
                      TRANSFORM = parent.frame()$TRANSFORM,TRANSFORM_X = parent.frame()$TRANSFORM_X, INVERSE_X = parent.frame()$INVERSE_X,
                      TRANSFORM_Y = parent.frame()$TRANSFORM_Y,INVERSE_Y = parent.frame()$INVERSE_Y,
                      FOD = parent.frame()$FOD,FOD_MODEL = parent.frame()$FOD_MODEL, FOD_SDRES_MIN = parent.frame()$FOD_SDRES_MIN, FOD_STDRES_MAX = parent.frame()$FOD_STDRES_MAX,
                      TRIMM = parent.frame()$TRIMM,
                      SOD = parent.frame()$SOD,SOD_MODEL = parent.frame()$SOD_MODEL, SOD_SDRES_MIN = parent.frame()$SOD_SDRES_MIN, SOD_STDRES_MAX = parent.frame()$SOD_STDRES_MAX,
                      LR_SD_RES_FACTOR = parent.frame()$LR_SD_RES_FACTOR, R2_MIN = parent.frame()$R2_MIN,
                      BATCH_HARMONIZATION = parent.frame()$BATCH_HARMONIZATION,
                      #CAL_CONC = parent.frame()$CAL_CONC,
                      GET_LR_STATUS = parent.frame()$GET_LR_STATUS,
                      nCORE = parent.frame()$nCORE,
                      GET_OUTPUT = parent.frame()$GET_OUTPUT, IMG_OUTPUT_DIR = parent.frame()$IMG_OUTPUT_DIR){

  data.table::setDT(dat)


  if(dim(dat)[1] < MIN_FEATURE) rlang::abort("'input_dat' has less rows than argument 'min_feature' required")
  if(dim(dat)[2] < 7) rlang::abort("'input_dat' needs minimum 7 columns defining the sample identification, the feature identification, the x and y values, the sample type of the data, the batch nr, the injection order and the dilution step.")

  # if (!"Batch" %in% names(COLNAMES)) {
  #   COLNAMES["Batch"] <- "Batch"
  #   dat[, Batch := "1"]
  # } else if("Batch" %in% names(COLNAMES) & is.na(COLNAMES[["Batch"]])){
  #   COLNAMES["Batch"] <- "Batch"
  #   dat[, Batch := "1"]
  # }


  if (!"Class" %in% names(COLNAMES)) {
    COLNAMES["Class"] <- "Class"
    dat[, Class := "1"]
  } else if("Class" %in% names(COLNAMES) & is.na(COLNAMES[["Class"]])){
    COLNAMES["Class"] <- "Class"
    dat[, Class := "1"]
  }

  if(!all(COLNAMES %in% colnames(dat)))rlang::abort("the provided names do not fit the colnames of the input data")


  data.table::set(dat, j = COLNAMES[["Sample_ID"]], value = as.character(dat[[COLNAMES[["Sample_ID"]]]]))
  data.table::set(dat, j = COLNAMES[["Feature_ID"]], value = as.character(dat[[COLNAMES[["Feature_ID"]]]]))
  data.table::set(dat, j = COLNAMES[["Batch"]], value = as.character(dat[[COLNAMES[["Batch"]]]]))
  data.table::set(dat, j = COLNAMES[["Injection_order"]], value = as.numeric(dat[[COLNAMES[["Injection_order"]]]]))
  data.table::set(dat, j = COLNAMES[["Class"]], value = as.character(dat[[COLNAMES[["Class"]]]]))
  data.table::set(dat, j = COLNAMES[["X"]], value = as.numeric(dat[[COLNAMES[["X"]]]]))
  data.table::set(dat, j = COLNAMES[["Y"]], value = as.numeric(dat[[COLNAMES[["Y"]]]]))
  data.table::set(dat, j = COLNAMES[["Sample_type"]], value = as.character(dat[[COLNAMES[["Sample_type"]]]]))


  data.table::setorderv(dat, c(COLNAMES[["Injection_order"]], COLNAMES[["Batch"]], COLNAMES[["Feature_ID"]], COLNAMES[["Sample_ID"]], COLNAMES[["Sample_type"]],COLNAMES[["Class"]], COLNAMES[["X"]], COLNAMES[["Y"]]))
  dat[ , IDintern := paste(dat[[COLNAMES[["Sample_ID"]]]], dat[[COLNAMES[["Feature_ID"]]]], dat[[COLNAMES[["Batch"]]]], dat[[COLNAMES[["Injection_order"]]]],sep = "_")]
  dat[ , groupIndices := .GRP ,by = c(COLNAMES[["Feature_ID"]], COLNAMES[["Batch"]])]
  dat <- dat[ , c( "IDintern", "groupIndices", COLNAMES[["Feature_ID"]],  COLNAMES[["Sample_ID"]], COLNAMES[["Sample_type"]],COLNAMES[["Class"]], COLNAMES[["Batch"]],COLNAMES[["Injection_order"]], COLNAMES[["X"]], COLNAMES[["Y"]]), with = F]

  # tests for Input parameter

  if(!(is.double(dat[[COLNAMES[["X"]]]]) & all(dat[[COLNAMES[["X"]]]] >= 0, na.rm = T))) rlang::abort("all values of 'column_X' and 'column_Y' need to be from type double and positive.")
  if(!(is.double(dat[[COLNAMES[["X"]]]]) & all(dat[[COLNAMES[["Y"]]]] >= 0, na.rm = T))) rlang::abort("all values of 'column_X' and 'column_Y' need to be from type double and positive.")
  if(!is.character(dat[[COLNAMES[["Feature_ID"]]]])) rlang::abort("'column_featureID' needs to be from type character")
  if(!is.character(dat[[COLNAMES[["Sample_ID"]]]])) rlang::abort("'column_sampleID' needs to be from type character")
  if(!is.numeric(dat[[COLNAMES[["Injection_order"]]]])) rlang::abort("'column_injectionOrder' needs to be from type numeric")
  if(!is.character(dat[[COLNAMES[["Batch"]]]])) rlang::abort("'column_Batch' needs to be from type character ")
  if(!is.character(dat[[COLNAMES[["Class"]]]])) rlang::abort("'column_sampleClass' needs to be from type character ")

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
  #if(!is.logical(CAL_CONC) | is.na(CAL_CONC)) rlang::abort("Argument 'calculate_concentration' needs to be from type logical")
  if(!is.logical(GET_LR_STATUS) | is.na(GET_LR_STATUS)) rlang::abort("Argument 'get_linearity_status_samples' needs to be from type logical")
  if(!nCORE > 0) rlang::abort("Argument 'nCore' needs to be positive")
  if(!MIN_FEATURE >= 3) rlang::abort("Argument 'min_feature' needs to be greater or equal than 3")
  if(!is.numeric(LR_SD_RES_FACTOR)) rlang::abort("Argument 'LR_sd_res_factor' needs to be from type integer and positive")
  if(!is.wholenumber(LR_SD_RES_FACTOR) | LR_SD_RES_FACTOR < 0) rlang::abort("Argument 'LR_sd_res_factor' needs to be from type integer and positive")

  if(!any(dat[[COLNAMES[["Sample_type"]]]] %in% CALIBRANTS)) rlang::abort("Argument 'sampleType_serial' was not found in column 'column_sampleType'")
  if(!data.table::between(x = R2_MIN, lower = 0, upper = 1)) rlang::abort("Argument 'R2_min' needs to be from type double and in the range between 0 and 1")
  if(!TYPE %in% c("untargeted", "targeted")) rlang::abort("Argument 'analysisType' need to be either 'untargeted' or 'targeted'")
  #if(TYPE %in% "untargeted")if(!(!is.null(DILUTION_FACTOR) | DILUTION_FACTOR != " " | !is.na(DILUTION_FACTOR))) rlang::abort("Argument 'dilution_factor ' must be provided for untargeted analysis")

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


  if(!all(is.na(SAMPLE)) & !SAMPLE %in% dat[[COLNAMES[["Sample_type"]]]]) rlang::abort("Argument 'sampleType_sample' was not found in column 'column_sampleType'")
  #if(!is.null(SAMPLE)) if(!any(colnames(dat) %in% Y_SAMPLE))rlang::abort("Column 'column_Y_sample' was not found in input data.")

  if(!all(is.na(QC))) if(!all( QC %in% dat[[COLNAMES[["Sample_type"]]]])) rlang::abort("Argument 'sampleType_QC' was not found in column 'column_sampleType'")
  #if(!is.numeric(DILUTION_FACTOR)) rlang::abort("Argument 'dilution_factor' needs to be from type numeric.")
  if(!is.logical(BATCH_HARMONIZATION) | is.na(BATCH_HARMONIZATION)) rlang::abort("Argument 'Batch_harmonization' needs to be from type logical")
  if(!is.null(NOISE)) if(!is.numeric(NOISE) | NOISE < 1) rlang::abort("Argument 'signal_blank_ration' needs to be from type integer and greater than 0.")
  if(!is.null(NOISE)) if(!BLANK %in% dat[[COLNAMES[["Sample_type"]]]]) rlang::abort("Argument 'sampleType_blank' was not found in column 'column_sample_type'")




  return(dat)
}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol



#' Transform and structure data for linearity modelling
#'
#' @description
#' This function prepares validated metabolomics data for downstream
#' linearity assessment by standardising column names, generating
#' internal identifiers, ordering observations, and optionally applying
#' mathematical transformations to concentration and intensity values.
#'
#' @details
#' The function assumes that the input data has already been validated
#' using \code{checkData()}.
#'
#' Processing steps include:
#' \itemize{
#'   \item Renaming of columns according to MSlineaR internal format
#'   \item Generation of dilution ordering index
#'   \item Creation of plotting metadata (e.g. point shape and colour)
#'   \item Feature-wise scaling of X values
#'   \item Optional transformation of X and Y variables
#' }
#'
#' Transformations (if enabled) are applied using user-defined functions
#' (e.g. log10, sqrt), and infinite values are automatically replaced with NA.
#'
#' @param dat A validated data.table returned by \code{checkData()}.
#'
#' @param TRANSFORM Logical. Whether transformations should be applied.
#'
#' @param TRANSFORM_X Character. Name of function used to transform X values.
#'
#' @param TRANSFORM_Y Character. Name of function used to transform Y values.
#'
#' @param COLNAMES Named list of column mappings used in the original dataset.
#'
#' @param TYPE Character. Type of analysis ("targeted" or "untargeted").
#'
#' @return A processed data.table containing:
#' \itemize{
#'   \item Standardised identifiers (IDintern, Feature_ID, Sample_ID)
#'   \item Processed concentration and intensity values (X, Y)
#'   \item Optional transformed variables (X_trans, Y_trans)
#'   \item Plotting metadata (DilutionPoint, pch, color)
#' }
#'
#' @examples
#' \dontrun{
#' dat_processed <- prepareData(dat_checked,
#'                              TRANSFORM = TRUE,
#'                              TRANSFORM_X = "log10",
#'                              TRANSFORM_Y = "log10")
#' }
#'
#' @keywords internal
#'
#' @import data.table
#' @importFrom rlang abort
prepareData <- function(dat,
                        TRANSFORM = parent.frame()$TRANSFORM,
                        TRANSFORM_X = parent.frame()$TRANSFORM_X,
                        TRANSFORM_Y = parent.frame()$TRANSFORM_Y,
                        #DILUTION_FACTOR = parent.frame()$DILUTION_FACTOR,
                        COLNAMES = parent.frame()$COLNAMES,
                        TYPE = parent.frame()$TYPE){

  if(dim(dat)[2] != 10) rlang::abort("data need to have 10 columns, please use funtion 'checkData' before to check all necessary input arguments")



  data.table::setDT(dat)
  processed <- data.table::copy(dat)
  # rename
  data.table::setnames(x = processed, old = colnames(processed), new = c( "IDintern","groupIndices", "Feature_ID","Sample_ID", "Sample.Type","Class", "Batch", "Injection_order",COLNAMES[["X"]], "Y"))
  data.table::set(processed, j = "Dilution", value = as.numeric(processed$Dilution))

  data.table::setorderv(processed, c(COLNAMES[["X"]],"Feature_ID", "Batch", COLNAMES[["X"]]))

  processed[ , ":="(Comment = NA, pch = data.table::fcase(!is.na(Y), 19), color = data.table::fcase(is.na(Y), "grey", default = "black"))]
  processed[ , ":="(#YNorm = Y / max(Y, na.rm = T) * 100,
    #Y_trans = log(Y),
    #X_trans = log(X),
    DilutionPoint = 1 : (.N),
    X = get(COLNAMES[["X"]]) * 3),
    #X = cumprod(c(1,get(COLNAMES[["X"]])[-1]/get(COLNAMES[["X"]])[- (.N)])) * 3),
    #groupIndices = .GRP),
    by = c("Feature_ID", "Batch")]

  #processed$X = processed[[COLNAMES[["X"]]]] * 100
  #processed$DilutionPoint <- processed$DilutionPoint + 1

  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_X)){
    processed$X_trans = get(TRANSFORM_X)(processed$X)
    processed$X_trans[is.infinite(processed$X_trans)] <- NA
  }

  if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)){
    processed$Y_trans = get(TRANSFORM_Y)(processed$Y)
    processed$Y_trans[is.infinite(processed$Y_trans)] <- NA
  }

  #processed <- processed[ , c("groupIndices", "IDintern","Feature_ID", "Batch", "X", "Y", "Comment", "YLog", "XLog", "DilutionPoint", "pch", "color"), with = F]

  return(processed)
}



#' Parallel execution helper for feature-wise analysis
#'
#' @description
#' This function provides a parallelised wrapper to apply a user-defined
#' function to groups of metabolomics features (defined by `groupIndices`).
#'
#' It is primarily used to accelerate computationally intensive steps
#' in the MSlineaR workflow, such as model fitting and linearity assessment
#' across large feature sets.
#'
#' @details
#' The function splits the input dataset into subsets based on unique
#' feature-group identifiers (`groupIndices`) and distributes these subsets
#' across multiple CPU cores using a PSOCK cluster.
#'
#' A progress bar is displayed during execution.
#'
#' Each subset is processed independently using the user-supplied function
#' `func`.
#'
#' @param nCORE Integer. Number of CPU cores to use for parallel processing.
#'
#' @param xs Integer vector. Indices of feature groups to be processed.
#' Each element corresponds to one unique value of `groupIndices`.
#'
#' @param func Function. User-defined function applied to each subset of data.
#' The function must accept a data.table as its first argument.
#'
#' @param inputData A data.frame or data.table containing at least the column
#' `groupIndices`, which defines feature-wise grouping.
#'
#' @param ... Additional arguments passed to `func`.
#'
#' @return
#' A list of results, where each element corresponds to the output of `func`
#' applied to one subset of `inputData` defined by `groupIndices`.
#'
#' @examples
#' \dontrun{
#' result <- my_fcn(
#'   nCORE = 4,
#'   xs = 1:10,
#'   func = function(dt) {
#'     mean(dt$Y, na.rm = TRUE)
#'   },
#'   inputData = dat
#' )
#' }
#'
#' @import data.table
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @keywords internal

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


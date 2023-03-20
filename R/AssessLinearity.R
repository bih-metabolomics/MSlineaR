#' Cleaning up a mass spectrometry data set according to the linear response of serial diluted/concentrated signals.
#'
#' @description `AssessLinearity()` cleans up a data set to only signals which
#' show a linear response. The function use two outlier detections and a partially
#' linear regression to identify the exact linear response range for each Signal.
#'
#'
#' @param input_data Long format data frame or data table combining the information
#' of samples, serial diluted/concentrated samples and optional repeated measured
#' QC samples for all batches.
#' @param column_sample_type String; Column name which distinguish between
#' samples, QC and serial samples.
#' @param sample_type_QC String; Identification of QC in `column_sample_type`.
#' @param sample_type_sample String;Identification of sample in `column_sample_type`.
#' @param sample_type_serial String;Identification of serial data in `column_sample_type`.
#' @param column_ID String; Column name, which identifies the measured Signals.
#' @param column_Batch String; Column name to identify the different batches.
#' Only necessary if there are more than 1 Batch.
#' @param column_X String; Column name of the independent variable,
#'  e.g. Concentration, Dilution,...
#' @param column_Y String; Column name of the dependent variable,
#'  e.g. Intensity, Area,..
#' @param transform Boolean; Should the data be transformed? Default is `TRUE`.
#' @param transform_x,transformy If `transform` is `TRUE`.
#' String; Which transformation should be used for the independent variable (concentration/dilution) or/and
#' the dependent variable of the serial diluted/concentrated samples?
#' Default for both is "log". If no transformation should be performed for one variable use NA.
#' @param first_outlier_detection,second_outlier_detection Boolean;
#' Should an outlier detection be performed before/after the exclusion of
#' nonlinear portions at the beginning/ end of the concentration/dilution range?
#' Default is `TRUE`
#' @param FOD_model,SOD_model Only necessary if `first_outlier_detection`/`second_outlier_detection` is `TRUE` (default).
#' A vector of statistical models, which should be used for the outlier detection.
#' Currently supported are linear regression ("linear"), quadratic regression ("quadratic")
#' and logistic regression ("logistic"). Default is a vector with all three models:
#' `FOD_model = c("linear", "logistic", "quadratic")`
#' @param FOD_sdres_min,SOD_sdres_min Only necessary if
#' `first_outlier_detection`/`second_outlier_detection` is `TRUE` (default).
#' Integer; Minimum residual standard deviation of a statistical model.
#' If the residual standard deviation is below this value, there will be no
#' outlier detection performed for this signal. Default to 1.
#' @param FOD_stdres_max,SOD_stdres_max Only necessary if
#' `first_outlier_detection`/`second_outlier_detection` is `TRUE` (default).
#' Integer; Maximum value for standardized residuals of a statistically model.
#' If a standardized residual is above this value.
#' this point will be considered as outlier and removed for further procedere. Default to 2.
#' @param R2_min
#' Numeric, ranges from 0 to 1; Minimum coefficient of determination,
#' which needs to be reached to consider the signal as linear. Default to 0.9
#' @param trimming Boolean; Should the data be trimmed? Default to `TRUE`.
#' @param min_feature Integer, ranging between 3 and maximum number of dilutions/concentrations.
#' Minimum number of points present in one serial diluted/concentrated series
#' marked as linear to consider this signal as linear. Default to 6, according to EMA guidelines2022.
#' @param LR_sd_res_factor Integer; points of serial diluted/concentrated series,
#' which are lower than `LR_sd_res_factor` times residual standard deviation are considered as linear.
#' Default to 2.
#' @param calculate_concentration Boolean; For targeted analysis. Should the
#' concentration of samples in regard to the linear regression equation be calculated? Default is `TRUE`.
#' @param get_linearity_status_samples Boolean; If `TRUE` (default) the samples
#' will be differentiated into linear and non linear according to their serial diluted/concentrated signal.
#' @param get_output Boolean; If `TRUE` (default), output files will be generated
#' and stored in output folder.
#' @param output_name Only necessary if `get_output` is `TRUE`. String; custom prefix
#' for the output files.
#' @param which_output Only necessary if `get_output` is `TRUE`. Vector of output files,
#'  which should be generated.Default all files will be generated and stored:
#' `which_output = c("R_object", "serial_list", "samples_all", "samples_filtered", "plot")`
#' @param output_dir Only necessary if `get_output` is `TRUE`.
#' custom path, where the output files should be stored.
#' @param nCORE Integer; Number of cores used for parallelization. Default to 1.
#' @param analysis_type String; was the analysis "targeted" or "untargeted"?
#' @param Batch_harmonization Boolean; If `TRUE` (default),
#' only serial diluted curves will be considered if they have a linear Range in
#' all Batches, otherwise they will be excluded. If set to `FALSE` all curves with
#' a linear range will be taken but flagged that in one or more batches linearity is missing
#'
#'
#' @importFrom  dplyr %>%
#' @return
#'
#'
#'
#' @export
#'
#' @examples
#'
AssessLinearity <- function(
    analysis_type = c("targeted", "untargeted", NULL)[3],
    #input_data
    input_data = NULL,
    column_sample_type = c("Sample.Type", "Type")[1],
    sample_type_QC = c("pooled QC", "QC")[1],
    sample_type_sample = "Sample",
    sample_type_serial = "Calibration Standard",
    column_ID = "Compound",
    column_Batch = "Batch",
    column_X = c("Concentration", "Dilution")[2],
    column_Y = c("Area","Height", "Intensity")[1],
    column_Y_sample = c("Area","Height", "Intensity")[1],

    #transform_data
    transform = c(TRUE, FALSE)[1],
    transform_x = "log",
    inverse_x = "exp",
    transform_y = "log",
    inverse_y = "exp",

    #first_outlier_detection
    first_outlier_detection = c(TRUE, FALSE)[1],
    FOD_model = c("logistic", "linear", "quadratic"),
    FOD_sdres_min = 1,
    FOD_stdres_max = 2,

    trimming = c(TRUE, FALSE)[1],

    #second_outlier_detection
    second_outlier_detection = c(TRUE, FALSE)[1],
    SOD_model = c("logistic", "linear", "quadratic"),
    SOD_sdres_min = 1,
    SOD_stdres_max = 2,

    #linear_range
    min_feature = 6,
    LR_sd_res_factor = 2,
    R2_min = 0.9,

    calculate_concentration = c(TRUE, FALSE)[1],

    get_linearity_status_samples = c(TRUE, FALSE)[1],

    Batch_harmonization = c(TRUE, FALSE)[1],

    nCORE = 1,

    get_output = c(TRUE, FALSE)[1],
    output_name = NULL,
    which_output = c("R_object", "allCurveSignal", "allCurveFeature",
                     "filteredCurveSignal", "allBioSampleSignal", "filteredBioSampleSignal", "plots")[1:7],
    output_dir = NULL
) {


  # define arguments
  TYPE = analysis_type
  DAT = input_data
  #QC =  sample_type_QC
  SAMPLE = sample_type_sample
  CALIBRANTS = sample_type_serial
  COLNAMES = c(ID = column_ID, Batch = column_Batch, Sample_type = column_sample_type, X = column_X,Y = column_Y)
  Y_SAMPLE = column_Y_sample
  TRANSFORM = transform
  TRANSFORM_X = transform_x
  INVERSE_X = inverse_x
  TRANSFORM_Y = transform_y
  INVERSE_Y = inverse_y
  FOD = first_outlier_detection
  FOD_MODEL = FOD_model
  FOD_SDRES_MIN = FOD_sdres_min
  FOD_STDRES_MAX = FOD_stdres_max
  #FOD_R2_MIN = FOD_R2_min
  TRIMM = trimming
  SOD = second_outlier_detection
  SOD_MODEL = SOD_model
  SOD_SDRES_MIN = SOD_sdres_min
  SOD_STDRES_MAX = SOD_stdres_max
  #SOD_R2_MIN = SOD_R2_min
  MIN_FEATURE = min_feature
  LR_SD_RES_FACTOR = LR_sd_res_factor
  R2_MIN = R2_min
  BATCH_HARMONIZATION = Batch_harmonization
  CAL_CONC = calculate_concentration
  GET_LR_STATUS = get_linearity_status_samples
  nCORE = nCORE
  GET_OUTPUT = get_output
  PREFIX = output_name
  OUTPUT_DIR = output_dir


  message("checking input arguments\n--------------------------------------------------------\n")
  dataOrigin <- checkData(input_data)


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


  # Creates full absolute paths
  if(GET_OUTPUT %in% TRUE){
    REPORT_OUTPUT_DIR <- R.utils::getAbsolutePath(OUTPUT_DIR)
    # check if output dir already exist or create it
    if (!dir.exists(REPORT_OUTPUT_DIR)){
      dir.create(REPORT_OUTPUT_DIR)
      print("Dir", REPORT_OUTPUT_DIR , "was created.")
    } else {
      print("Dir already exists!")
    }
  }


  if(TRANSFORM %in% TRUE){
    X = ifelse(is.na(TRANSFORM_X), "X", "X_trans")
    Y = ifelse(is.na(TRANSFORM_Y), "Y", "Y_trans")

  } else{
    X = "X"
    Y = "Y"
  }

  Xraw = X
  Yraw = Y

  #wording

  if(TYPE %in% "targeted"){
    Compounds = "Compound(s)"
    Dilutions = "Concentrations"
    Series = "Concentration Curves"
    Signals = "Signals"
  } else{
    Compounds = "molecular feature"
    Dilutions = "Dilutions"
    Series = "QC dilution Curves"
    Signals = "Signals"
  }

  # if (!"REPLICATE" %in% names(COLNAMES)) {
  #   COLNAMES["REPLICATE"] <- "REPLICATE"
  #   dataOrigin[, REPLICATE := "1"]
  # }



  processingFeatureCal <- dataOrigin[get(COLNAMES[["Sample_type"]]) %in% CALIBRANTS]
  processingFeature <- data.table::copy(processingFeatureCal)


  #### normalizing, centralizing, log transforming ####
  message("preparing serial diluted QC data\n--------------------------------------------------------\n")

  dataPrep <- prepareData(processingFeature)


  processingFeature <- data.table::copy(dataPrep)

  step = 1
  countList <- countMinimumValue(DAT = processingFeature,step = step,y = Y)
  processingFeature <- data.table::data.table(countList[[1]])
  processingGroup <- data.table::data.table(countList[[2]])


  # save key numbers
  nCompounds <- data.table::uniqueN(processingGroup, by = c("ID"))
  nReplicates <- data.table::uniqueN(processingFeature, by = c("Batch"))
  nDilutions <- data.table::uniqueN(processingFeature[, 'DilutionPoint' := 1:.N, by = c("ID", "Batch")], by = c("DilutionPoint"))
  nSeries <- data.table::uniqueN(processingGroup, by = c("groupIndices"))
  nPeaks <- data.table::uniqueN(processingFeature, by = c("IDintern"))

  message(
    "Data set with ", nCompounds, " ",Compounds, ", ",
    nReplicates, " Batch(es) and ",
    nDilutions, " ", Dilutions," -> ",
    format.default(nSeries, big.mark = ",", scientific = F), " ",Series,".",
    "\n--------------------------------------------------------\n"
  )

  # check length of points

  nCompoundsNew <- data.table::uniqueN(processingGroup[enoughPeaks_1 %in% TRUE], by = c("ID"))
  nReplicatesNew <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("Batch"))
  nSeriesNew <- data.table::uniqueN(processingGroup[enoughPeaks_1 %in% TRUE], by = c("groupIndices"))

  message(paste0("Removed ", nSeries - nSeriesNew, " ", Series, " with less than ", MIN_FEATURE, " Signals.","\n",
                 "Remaining Data set with ", nCompoundsNew, " ", Compounds," and ", nReplicatesNew, " Batch(es) -> ", nSeriesNew, " ", Series),"\n--------------------------------------------------------\n")

  stopifnot(exprs = {
    "all Compounds were removed" = nCompoundsNew - data.table::uniqueN(processingGroup[enoughPeaks_1 %in% FALSE], by = c("ID")) > 0
    "all Dilution/Concentration-Series were removed" = nSeriesNew - data.table::uniqueN(processingGroup[enoughPeaks_1 %in% FALSE], by = c("groupIndices")) > 0
  })

  cutoff <- processingFeature[groupIndices %in% processingGroup[enoughPeaks_1 %in% FALSE, groupIndices]]
  processingFeature <- processingFeature[groupIndices %in% processingGroup[enoughPeaks_1 %in% TRUE, groupIndices]]

  step = step + 1

  # assert_that(nrow(processList$processing) == rawCompounds*n_dilution*nReplications)

  #### first outlier detection ####
  if(FOD %in% TRUE){
    message("First Outlier Detection\n")

    # prints recorded time
    startTime = Sys.time()

    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    #on.exit(parallel::stopCluster(cl))
    #options(future.globals.onReference = "error")
    dataFOD <- my_fcn(
      #cl,
      #exportObjects = c("chooseModel", "X","Y", "abbr", "R2min"),
      xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      x = X,
      y = Y,
      func = chooseModel,
      abbr = "FOD",
      #R2_MIN = FOD_R2_MIN,
      model = FOD_MODEL,
      SDRES_MIN = FOD_SDRES_MIN,
      STDRES = FOD_STDRES_MAX
    ) |> unlist(recursive = F)

    #parallel::stopCluster(cl)
    closeAllConnections()

    message("FOD: ",format(round(difftime(Sys.time(), startTime),2)))
    Y = "Y_FOD"

    dataFODModel <- tibble::tibble(
      groupIndices = as.integer(names(purrr::map(dataFOD, 1))),
      ModelName = purrr::map(dataFOD, 1) %>% unlist(use.names = F),
      Model = purrr::map(dataFOD, 2),
      #RMSE = purrr::map(dataFOD,3) %>% unlist(use.names = T),
      #R2 = purrr::map(dataFOD, 3) %>% unlist(use.names = F)
    )

    dataFOD = purrr::map(dataFOD,3)|> plyr::ldply(.id = NULL)

    processingFeature <- data.table::data.table(dplyr::full_join(dataFOD, processingFeature[!IDintern %in% dataFOD$IDintern], by = colnames(processingFeature)))
    processingGroup <- dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,'OutlierFOD' :=any(OutlierFOD %in% TRUE), groupIndices][,.(groupIndices, OutlierFOD)]), by = c("groupIndices"))

    message("An Outlier were found for ", data.table::uniqueN(processingFeature |> dplyr::filter(OutlierFOD %in% TRUE) %>% dplyr::select(groupIndices)), " ",Series,".\n")

    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- countList[[1]]
    processingGroup <- full_join(processingGroup, countList[[2]], by = c("groupIndices", "ID", "Batch"))

    # check length of points
    checkData <- checkLength()
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1


  }
  #### trim ####
  if(TRIMM %in% TRUE){

    message("Trim data: first Dilution should have the smallest Intensity and last point should have the biggest.")

    startTime = Sys.time()
    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    dataTrim <- my_fcn(
     # cl,
      xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      func = trimEnds,
      x = X,
      y = Y
    ) |>
      plyr::ldply(.id = NULL)

    #parallel::stopCluster(cl)
    closeAllConnections()

    message("Trimming: ",format(round(difftime(Sys.time(), startTime),2)))

    Y = "Y_trim"

    processingFeature <- data.table::data.table(dplyr::full_join(dataTrim,
                                                                 processingFeature[!IDintern %in% dataTrim$IDintern],
                                                                 by = colnames(processingFeature)))
    processingGroup <- dplyr::full_join(processingGroup,
                                        unique(data.table::copy(processingFeature)[,'trim' :=any(trim %in% TRUE),groupIndices][,.(groupIndices, trim)]),
                                        by = "groupIndices")

    TrimGroups <- data.table::uniqueN(dataTrim |> dplyr::filter(trim %in% TRUE) %>% dplyr::select(groupIndices))
    TrimFeatures <- data.table::uniqueN(dataTrim |> dplyr::filter(trim %in% TRUE) |> dplyr::select(IDintern))


    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- dplyr::full_join(countList[[1]], processingFeature[!IDintern %in% countList[[1]]$IDintern], by = colnames(processingFeature))
    processingGroup <- full_join(processingGroup, countList[[2]], by = c("groupIndices", "ID", "Batch"))

    message("In total ", TrimFeatures," ",  Signals," in ", TrimGroups, " ",Series," were trimmed.\n")

    # check length of points
    checkData <- checkLength()
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1
    # assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))

  }

  #### second outlier detection ####

  if(SOD %in% TRUE){

    message("Second Outlier Detection\n")


    processingFeature <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE, groupIndices]]

    startTime = Sys.time()
    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    dataSOD <- my_fcn(
      #cl,
      xs = 1:data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      x = X,
      y = Y,
      func = chooseModel,
      abbr = "SOD",
      #R2_MIN = SOD_R2_MIN,
      model = SOD_MODEL,
      SDRES_MIN = SOD_SDRES_MIN,
      STDRES = SOD_STDRES_MAX
    ) |> unlist(recursive = F)

    closeAllConnections()

    #rm(cl)


    endTime = Sys.time()
    message("SOD: ",format(round(difftime(Sys.time(), startTime),2)))

    Y = "Y_SOD"


    dataSODModel <- tibble::tibble(
      groupIndices = as.integer(names(purrr::map(dataSOD, 1))),
      ModelName = purrr::map(dataSOD, 1) %>% unlist(use.names = F),
      Model = purrr::map(dataSOD, 2),
      #RMSE = purrr::map(dataFOD,3) %>% unlist(use.names = T),
      #R2 = purrr::map(dataSOD, 3) %>% unlist(use.names = F)
    )

    #data.table::setDT(dataSODModel)
    dataSOD = purrr::map(dataSOD,3)|> plyr::ldply(.id = NULL)

    processingFeature <- data.table::data.table(dplyr::full_join(dataSOD, processingFeature[!IDintern %in% dataSOD$IDintern], by = colnames(processingFeature)))
    processingGroup <- dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,'OutlierSOD' :=any(OutlierSOD %in% TRUE), groupIndices][,.(groupIndices, OutlierSOD)]), by = c("groupIndices"))

    message("An Outlier were found for ", data.table::uniqueN(processingFeature |> dplyr::filter(OutlierSOD %in% TRUE) %>% dplyr::select(groupIndices)), " ",Series,".\n")

    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- countList[[1]]
    processingGroup <- full_join(processingGroup, countList[[2]], by = c("groupIndices", "ID", "Batch"))

    # check length of points
    checkData <- checkLength()
    processingFeature <- checkData[[1]]
    cutoffnew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1

  }

  #### positive association ####
{
    message("Check if data curves have a positive slope")



    startTime = Sys.time()
    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    dataTrimPos <- my_fcn(
      #cl,
      xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      func = trim_pos_associated,
      x = X,
      y = Y,
      MIN_Feature = MIN_FEATURE
    ) |>
      plyr::ldply(.id = NULL)

    closeAllConnections()

    message("TrimmPos: ",format(round(difftime(Sys.time(), startTime),2)))

    Y = "Y_trimPos"

    processingFeature <- data.table::data.table(dplyr::full_join(dataTrimPos,
                                                                 processingFeature[!IDintern %in% dataTrimPos$IDintern],
                                                                 by = colnames(processingFeature)))
    processingGroup <- dplyr::full_join(processingGroup,
                                        unique(data.table::copy(processingFeature)[,'trimPos' :=any(trimPos %in% TRUE),groupIndices][,.(groupIndices, trimPos)]),
                                        by = "groupIndices")

    TrimPosGroups <- data.table::uniqueN(dataTrimPos |> dplyr::filter(trimPos %in% TRUE) %>% dplyr::select(groupIndices))
    TrimPosFeatures <- data.table::uniqueN(dataTrimPos |> dplyr::filter(trimPos %in% TRUE) |> dplyr::select(IDintern))


    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- dplyr::full_join(countList[[1]], processingFeature[!IDintern %in% countList[[1]]$IDintern], by = colnames(processingFeature))
    processingGroup <- full_join(processingGroup, countList[[2]], by = c("groupIndices", "ID", "Batch"))

    message("In total ", TrimPosFeatures," ",  Signals," in ", TrimPosGroups, " ",Series," were removed.\n")

    # check length of points
    checkData <- checkLength()
    processingFeature <- checkData[[1]]
    cutoffnew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1
    # assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))



}




  ##### find linear Range ####
{  message("Determining linear range\n--------------------------------------------------------\n")

  processingFeature <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE, groupIndices]]

  #dataLin <-  data.table::copy(processingFeature)[groupIndices %in% dataSODModel[aboveMinCor %in% TRUE, groupIndices]]
  #dataLin <-  data.table::copy(processingFeature)[groupIndices %in% dataSODModel[ , groupIndices]]
  #dataLin$fittingModel <- sapply(dataLin$groupIndices, function(i) dataSODModel$Model[dataSODModel$groupIndices %in% i])

  startTime = Sys.time()
  #cl <- parallel::makeCluster(getOption("cl.cores", nCORE), type = "PSOCK")


  dataLinearRange <- my_fcn(
    #cl = myCluster,
    xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    func = findLinearRange,
    x = X,
    y = Y,
    min_feature = MIN_FEATURE,
    sd_res_factor = LR_SD_RES_FACTOR
  )

  closeAllConnections()
  message("FindLinear Range: ",format(round(difftime(Sys.time(), startTime),2)))

  Y = "Y_LR"

  dataLRFeature = data.table::data.table(purrr::map(dataLinearRange,1)|> plyr::ldply(.id = NULL))
  dataLRGroup = data.table::data.table(purrr::map(dataLinearRange,2)|> plyr::ldply(.id = NULL))
  dataLRGroup$aboveR2 = dataLRGroup$R2 >= R2_MIN
  dataLRGroup$LRFlag[dataLRGroup$aboveR2 %in% FALSE] = "small R2"
  dataLRGroup$enoughPointsWithinLR[dataLRGroup$aboveR2 %in% FALSE] = FALSE


  dataLRFeaturesub <- dataLRFeature[groupIndices %in% dataLRGroup[aboveR2 %in% FALSE, groupIndices] & color %in% "darkseagreen", .(IDintern, color, IsLinear, Comment)]
  dataLRFeature[dataLRFeaturesub,':=' ( color = "black",
                                        IsLinear = FALSE,
                                        Comment = paste0(Comment, "_small R2")) , on = "IDintern"]


  processingFeature <- data.table::data.table(dplyr::full_join(dataLRFeature, processingFeature[!IDintern %in% dataLRFeature$IDintern], by = colnames(processingFeature)))
  processingGroup <-  dplyr::full_join(processingGroup, dataLRGroup, by = "groupIndices")

  message("For ", data.table::uniqueN(processingGroup |> dplyr::filter(aboveR2 %in% TRUE) %>% dplyr::select(groupIndices)), " ",Series,
          " a linear Range with a minimum of ", MIN_FEATURE, " Points and an R^2 higher or equal ", R2_MIN," were found\n")

  countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
  processingFeature <- countList[[1]]
  processingGroup <- full_join(processingGroup, countList[[2]], by = c("groupIndices", "ID", "Batch"))



  # check length of points
  checkData <- checkLength()
  processingFeature <- checkData[[1]]
  cutoffNew <- checkData [[2]]

  cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))
}

######
  message("check batch reproducibility")

  {
    data.table::setDT(processingFeature)
    data.table::setDT(processingGroup)

    processingGroup <-   processingGroup[, enoughPeak_allBatches := all(get(paste0("enoughPeaks_", step)) %in% TRUE), by = ID]
    conflictBatches <- processingGroup[enoughPeak_allBatches %in% FALSE]
    processingFeature$Y_allBatches <- processingFeature$Y_LR

    if(BATCH_HARMONIZATION == TRUE){
      processingFeature[groupIndices %in% conflictBatches$groupIndices, Y_allBatches := NA]
      processingFeature[groupIndices %in% conflictBatches$groupIndices, color := "grey"]
      }

    processingFeature[groupIndices %in% conflictBatches$groupIndices,
                      Comment := paste(processingFeature[groupIndices %in% conflictBatches$groupIndices, Comment],"_notEnoughPeaksInAllBatches")]

    message(data.table::uniqueN(conflictBatches$ID), " ", Compounds, " were excluded, because they do not show repeatable ", Series, ".\n")
    message("Final Data set with ", data.table::uniqueN(processingGroup[enoughPeak_allBatches %in% TRUE, ID]), " ", Compounds,
            " and ", data.table::uniqueN(processingGroup[enoughPeak_allBatches %in% TRUE, Batch]), " Batch(es) -> ",
            data.table::uniqueN(processingGroup[enoughPeak_allBatches %in% TRUE, groupIndices]), " ", Series,"\n--------------------------------------------------------\n")

  }

  ### Back calculation Concentration
  if(TYPE %in% "targeted"){
    message("Back calculation of concentration for the concentration series signals")
    Y = ifelse(TRANSFORM %in% TRUE & !is.na(TRANSFORM_Y), "Y_trans", COLNAMES[["Y"]])
    processingFeature = getConc(dats = processingFeature, datCal = processingGroup, Y, INVERSE_Y)

  }
  processingFeature <- processingFeature[!is.na(IDintern)]


  #### biological samples ####
  if(CAL_CONC %in% TRUE | GET_LR_STATUS %in% TRUE){

  message("prepare biological samples")
    SampleFeature <- dataOrigin[get(COLNAMES[["Sample_type"]]) %in% SAMPLE]

    if(TRANSFORM %in% TRUE & !is.na(TRANSFORM_Y)){
      SampleFeature$Y_trans = get(TRANSFORM_Y)(SampleFeature[[Y_SAMPLE]])
      SampleFeature$Y_trans[is.infinite(SampleFeature$Y_trans)] <- NA
      Y_SAMPLE = "Y_trans"
    }

    if(TYPE %in% "targeted" & CAL_CONC %in% TRUE){
      message("Calculate Concentration for biological Samples")
    SampleFeature <- getConc(dats = SampleFeature, datCal = processingGroup, y = Y_SAMPLE, INVERSE_Y)
    }

    if(GET_LR_STATUS %in% TRUE){

      SampleFeature  <- getLRstatus(dats = SampleFeature, datCal = processingGroup,y =  column_Y_sample)
        #data_Sample$Status_LR <-  data.table::between(lower = data_Sample$LRStartY, x = data_Sample[, c("Area.CorrDrift")], upper = data_Sample$LREndY, NAbounds = NA)

    }
  }




  # assert_that(n_distinct(processList$processing$groupIndices) == rawCompounds)
  message("prepare output files\n--------------------------------------------------------\n")

  cutoff <- dplyr::full_join(cutoff,
                             processingFeature,
                             by = colnames(cutoff))


  processingFeature <- dplyr::full_join(processingFeature, cutoff, by = colnames(cutoff))

  assertthat::are_equal(nrow(processingFeature), nrow(processingFeatureCal))

  #output:
  data.table::setnames(processingFeature, c("ID","Sample.Type", "Batch", "Y", "X"), c(COLNAMES[["ID"]],COLNAMES[["Sample_type"]], COLNAMES[["Batch"]], COLNAMES[["Y"]], COLNAMES[["X"]] ))
  data.table::setnames(processingGroup, c("ID", "Batch"), c(COLNAMES[["ID"]], COLNAMES[["Batch"]]))


  #1)full table dilution/concentration curves - Signal based
  output1 <- processingFeature

  #2) full table dilution/concentration curves - Feature based
  output2 <- processingGroup

  #3) filtered table dilution/concentration curves - Signal based (high quality)
  output3 <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_",step)) %in% TRUE, groupIndices]]
  output3 <- output3 |> dplyr::group_by(groupIndices) |> dplyr::filter(!all(is.na(Y_allBatches)))

  #4) full table biological Samples - Signal based
  output4 <- SampleFeature

  #5) filtered table biological Samples - Signal based (high quality)
  output5 <- SampleFeature[Status_LR %in% TRUE]

  #6) barplot summary per dilution/concentration



  #7) scatter plot



  #8) barplot summary for biological samples






  data.table::setnames(processingGroup, c("ID", "Batch", "Y", "X"), c(COLNAMES[["ID"]],COLNAMES[["Batch"]], COLNAMES[["Y"]], COLNAMES[["X"]] ))

  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  processList <- list(
    "fullTableCurvesSignal" = output1,
    "fullTableCurvesFeature" = output2,
    "filteredTableCurvesSignal" = output3,
    "fullTableBioSamplesSignal" = output4,
    "filteredTableBioSamplesSignal" = output5,
    "dataModel_FOD" = dataFODModel,
    "dataModel_SOD" = dataSODModel,
    "Parameters" = list(
      data = DAT,
      COLNAMES = COLNAMES,
      nCORE = nCORE,
      LOG_TRANSFORM = LOG_TRANSFORM,
      R2min = R2min,
      res = res,
      #calcFor = calcFor,
      #SRES = SRES,
      #R2FOD = R2FOD,
      #R2SOD = R2SOD,
      MIN_FEATURE = MIN_FEATURE,
      X = X,
      Y = Y
    )
  )

  if(GET_OUTPUT %in% TRUE){
    if("allCurveSignal" %in% which_output) write.csv(output1, file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "allCurveSignal.csv" , sep = "_")))
    if("allCurveFeature" %in% which_output) write.csv(output2, file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "allCurveFeature.csv" , sep = "_")))
    if("filteredCurveSignal" %in% which_output) write.csv(output3, file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "filteredCurveSignal.csv" , sep = "_")))
    if("allBioSampleSignal" %in% which_output) write.csv(output4, file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "allBioSamplesSignal.csv" , sep = "_")))
    if("filteredBioSampleSignal" %in% which_output) write.csv(output5, file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "filteredBioSamplesSignal.csv" , sep = "_")))

  }


  # "dataOrigin_1" = dataOrigin,
  #  "dataPrep_2" = dataPrep,
  # "dataFOD_3" <- dataFOD,

  #"dataTrim_4" = dataTrim,
  #  "dataSOD_5" = dataSOD,

   # "dataLinearRange_7" = dataLinearRange,
    #"dataLR_8" = dataLRFeature,
    #"summaryFFDS" = processing,
    #"summaryFDS" = processingGroup,




  #if (!empty(dataFOD)) processList$"dataFOD_3" <- dataFOD
  #if (!is.null(dataSOD) | !empty(dataSOD)) processList$"dataSOD_7" <- dataSOD
  #if (!is.null(dataConsSOD) | !empty(dataConsSOD)) processList$"dataConsSOD_8" <- dataConsSOD
  #if (!is.null(dataModelSOD) | !empty(dataModelSOD)) processList$"dataModelSOD_9" <- dataModelSOD


  return(processList)
}

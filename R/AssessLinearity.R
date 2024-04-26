#' Cleaning up a mass spectrometry data set according to the linear response of serial diluted/concentrated signals.
#'
#' @description `AssessLinearity()` cleans up a data set of serial concentrated or diluted Compounds to only signals which
#' show a linear response. The function use two outlier detection steps and a partially
#' linear regression to identify the exact linear response range for each Compound.
#'
#' @param inputData_feature Long format data frame or data table combining the information
#' of biological samples, serial diluted/concentrated samples and optional repeated measured
#' QC samples for all batches.
#'  @param inputData_sample Long format data frame or data table combining the information
#' of biological samples, serial diluted/concentrated samples and optional repeated measured
#' QC samples for all batches.
#' @param column_sampleType String; Column name which distinguish between
#' samples, QC and serial samples.
#' @param sampleType_QC String; Identification of QC in `column_sampleType`.
#' @param sampleType_sample String;Identification of biological samples in `column_sampleType`.
#' @param sampleType_serial String;Identification of serial data in `column_sampleType`.
#' @param sampleType_blank String;Identification of blank data in `column_sampleType`.
#' @param column_sampleID String; Column name, which uniquely identifies the measured Signals.
#' @param column_batch String; Column name to identify the different batches.
#' Only necessary if there are more than 1 Batch.
#' @param column_X String; Column name of the independent variable,
#'  e.g. Concentration, Dilution,... . Highest number should be the highest concentration or lowest dilution.
#'  For untargeted data use a vector with consecutive numbers, e.g. c(1, 2, 3 , ...) with 1 beeing the highest dilution
#' @param column_Y String; Column name of the dependent variable,
#'  e.g. Intensity, Area,..
#' @param column_sampleClass String, Column name used for statistical classes
#' @param transform Boolean; Should the data be transformed? Default is `TRUE`.
#' @param transform_x,transform_y If `transform` is `TRUE`.
#' String; Which transformation should be used for the independent variable (concentration/dilution) or/and
#' the dependent variable of the serial diluted/concentrated samples?
#' Default for both is "log". If no transformation should be performed for one variable use NULL.
#' @param signal_blank_ratio Numeric,
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
#' which are less than `LR_sd_res_factor` times residual standard deviation are considered as linear.
#' Default to 3.
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
#' @param analysisType String; was the analysis "targeted" or "untargeted"?
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
MS_AssessLinearity <- function(
    analysisType = c("targeted", "untargeted", NULL)[1],

    #input_data
    inputData_feature = NULL,
    inputData_sample = NULL,
    column_sampleType = c("Sample.Type", "Type")[1],
    sampleType_QC = c("QCpool", "Reference QC", NULL)[1],
    sampleType_sample = "Sample",
    sampleType_serial = "Dilution",
    sampleType_blank = c("Blank", NULL)[1],
    signal_blank_ratio = 5,
    column_sampleID = "Sample.Identification",
    column_featureID = "Compound",
    column_injectionOrder = "Sequence.Position",
    column_batch = "Batch",
    column_X = c("Concentration", "Dilution")[2],
    column_Y = c("Area","Height", "Intensity")[1],
    #column_Y_sample = c("Area","Height", "Intensity")[1],
    column_sampleClass = c("Group", NULL)[1],
    #dilution_factor = NULL,

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
    LR_sd_res_factor = 3,
    R2_min = 0.9,

    Batch_harmonization = c(TRUE, FALSE)[1],

    #biological signals
    #calculate_concentration = c(TRUE, FALSE)[1],
    get_linearity_status_samples = c(TRUE, FALSE)[1],

    #filtering



    nCORE = 1,

    #output
    get_output = c(TRUE, FALSE)[1],
    output_name = NULL,
    which_output = c("all","R_object", "DilutionCurves", "BiologicalSamples", "Plots")[1],
    output_dir = NULL
) {


  # define arguments
  TYPE = analysisType
  DAT_F = inputData_feature
  DAT_S = inputData_sample
  QC =  sampleType_QC
  BLANK = sampleType_blank
  SAMPLE = sampleType_sample
  SAMPLE_ID = column_sampleID
  FEATURE_ID = column_featureID
  CALIBRANTS = sampleType_serial
  COLNAMES = c(Feature_ID = column_featureID,
               Sample_ID = column_sampleID,
               Injection_order = column_injectionOrder,
               Batch = column_batch,
               Sample_type = column_sampleType,
               X = column_X,
               Y = column_Y,
               Class = column_sampleClass)
  NOISE = signal_blank_ratio
  #Y_SAMPLE = column_Y_sample
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
  #CAL_CONC = calculate_concentration
  GET_LR_STATUS = get_linearity_status_samples
  nCORE = nCORE
  GET_OUTPUT = get_output
  PREFIX = output_name
  OUTPUT_DIR = output_dir

  #my_log <- file(paste0(OUTPUT_DIR,"/my_log.txt"), open = "a")
  #sink(my_log, append = TRUE, type = "output", split = TRUE)
  #sink(file = my_log, type = "message", append = TRUE)
  # sink.number()
  # #you can close all open sinks with
  # while (sink.number()>0) sink()


  # Create temp file location
  tmp <- file.path(OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_MSlineaR.log"))

  # Open log
  if (!dir.exists(OUTPUT_DIR)){
    rlang::inform(paste("Dir", OUTPUT_DIR , "was created."))
  } else {
    rlang::inform(paste("Dir", OUTPUT_DIR,"already exists!"))
  }
  lf <- logr::log_open(tmp, autolog = TRUE, show_notes = FALSE)

  # Send message to log
  #log_print(

  logr::sep("checking input arguments")
  # rlang::inform("
  #               --------------------------------------------------------
  #               \tchecking input arguments
  #               --------------------------------------------------------\n")

  #)



  dataOrigin_F <- data.table::copy(DAT_F)
  dataOrigin_S <- data.table::copy(DAT_S)

  #combine input tables
  testthat::expect_setequal(unique(dataOrigin_F[[SAMPLE_ID]]), unique(dataOrigin_S[[SAMPLE_ID]]))
  testthat::expect_setequal(unique(dataOrigin_F[[COLNAMES[["Batch"]]]]), unique(dataOrigin_S[[COLNAMES[["Batch"]]]]))

dataOrigin <- dplyr::full_join(dataOrigin_F, dataOrigin_S, by = intersect(colnames(dataOrigin_S), colnames(dataOrigin_F)))




dataReduced <- data.table::copy(dataOrigin)

  dataReduced <- checkData(dat = dataReduced) # function in prep.R


  # progressbar
  pbapply::pboptions(type = "timer", char = "[=-]", style = 5)

  progressr::handlers(global = NA)
  progressr::handlers(on_missing = "ignore")
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
    IMG_OUTPUT_DIR <- R.utils::getAbsolutePath(file.path(REPORT_OUTPUT_DIR, "img"))
    # check if output dir already exist or create it
    if (!dir.exists(REPORT_OUTPUT_DIR)){
      dir.create(REPORT_OUTPUT_DIR)
      logr::put(paste("Dir", REPORT_OUTPUT_DIR , "was created."))
      dir.create(file.path(IMG_OUTPUT_DIR))

    } else {
      #logr::put(paste("Dir", REPORT_OUTPUT_DIR,"already exists!"))
      if (!dir.exists(IMG_OUTPUT_DIR)){
        dir.create(IMG_OUTPUT_DIR)
      }}


      #check if pdfs are still open
      imgfileName = c("Summary_Barplot_QC", "Summary_Barplot_Samples", "Summary_Barplot_All", "Calibrationplot")

      if (any(c(file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[1],".pdf"))),
               file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[2],".pdf"))),
               file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[3],".pdf"))),
               file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[4],".pdf")))))){


        img.exist <- which(c(file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[1],".pdf"))),
                            file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[2],".pdf"))),
                            file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[3],".pdf"))),
                            file.exists(file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[4],".pdf")))))


        check <- sapply(img.exist, function(i) {file.rename(from = file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[i],".pdf")),
                                                              to = file.path(IMG_OUTPUT_DIR,paste0(Sys.Date(),"_", imgfileName[i],".pdf")) )
        })

        if(!all(check)) rlang::abort("Please close all pdfs first.")
      }








    }
  #}

Xorigin <- "X"
Yorigin <- "Y"
  if(TRANSFORM %in% TRUE){
    X <- ifelse(is.null(TRANSFORM_X), "X", "X_trans")
    Y <- ifelse(is.null(TRANSFORM_Y), "Y", "Y_trans")

  } else{
    X <- "X"
    Y <- "Y"
  }

  Xraw <- X
  Yraw <- Y

  #wording

  if(TYPE %in% "targeted"){
    Compounds <- "Compound(s)"
    Dilutions <- "Concentrations"
    Series <- "Concentration Curves"
    Signals <- "Signals"
  } else{
    Compounds <- "molecular feature"
    Dilutions <- "Dilutions"
    Series <- "QC dilution Curves"
    Signals <- "Signals"
  }

  processingFeatureCal <- dataReduced[get(COLNAMES[["Sample_type"]]) %in% CALIBRANTS]
  processingFeature <- data.table::copy(processingFeatureCal)


  #### normalizing, centralizing, log transforming ####
  logr::sep("preparing serial diluted QC data")
 # rlang::inform("
 #               --------------------------------------------------------
 #                \tpreparing serial diluted QC data
 #               --------------------------------------------------------\n")
# function in prep.R
  dataPrep <- prepareData(processingFeature, TRANSFORM, TRANSFORM_X, TRANSFORM_Y, COLNAMES, TYPE)

  processingFeature <- data.table::copy(dataPrep)

  step <- 1
  countList <- countMinimumValue(DAT = processingFeature,MIN_FEATURE,step = step,y = Y)
  processingFeature <- data.table::data.table(countList[[1]])
  processingGroup <- data.table::data.table(countList[[2]])


  # save key numbers
  nCompounds <- data.table::uniqueN(processingGroup, by = c("Feature_ID"))
  nReplicates <- data.table::uniqueN(processingFeature, by = c("Batch"))
  nDilutions <- data.table::uniqueN(processingFeature, by = c("DilutionPoint")  )
  nSeries <- data.table::uniqueN(processingGroup, by = c("groupIndices"))
  nPeaks <- data.table::uniqueN(processingFeature, by = c("IDintern"))

  logr::put(paste0(
    "Data set with ", nCompounds, " ",Compounds, ", ",
    nReplicates, " Batch(es) and ",
    nDilutions, " ", Dilutions," -> ",
    format.default(nSeries, big.mark = ",", scientific = F), " ",Series,"."
  ))

  # check length of points

  nCompoundsNew <- data.table::uniqueN(processingGroup[enoughPeaks_1 %in% TRUE], by = c("Feature_ID"))
  nReplicatesNew <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("Batch"))
  nSeriesNew <- data.table::uniqueN(processingGroup[enoughPeaks_1 %in% TRUE], by = c("groupIndices"))

  logr::put(paste0("Removed ", nSeries - nSeriesNew, " ", Series, " with less than ", MIN_FEATURE, " Signals.",
                 "Remaining Data set with ", nCompoundsNew, " ", Compounds," and ", nReplicatesNew, " Batch(es) -> ", nSeriesNew, " ", Series,"."))

  stopifnot(exprs = {
    "all Compounds were removed" = nCompounds - data.table::uniqueN(processingGroup[enoughPeaks_1 %in% FALSE], by = c("Feature_ID")) > 0
    "all Dilution/Concentration-Series were removed" = nSeries - data.table::uniqueN(processingGroup[enoughPeaks_1 %in% FALSE], by = c("groupIndices")) > 0
  })

  cutoff <- processingFeature[groupIndices %in% processingGroup[enoughPeaks_1 %in% FALSE, groupIndices]]
  processingFeature <- processingFeature[groupIndices %in% processingGroup[enoughPeaks_1 %in% TRUE, groupIndices]]

  step <- step + 1

  # assert_that(nrow(processList$processing) == rawCompounds*n_dilution*nReplications)

  #### remove Features < signal to blank ration ####
  if(!is.null(NOISE)){

    logr::sep("remove Features below signal to blank ratio")

    # rlang::inform("
    #            --------------------------------------------------------
    #             \tremove Features below signal to blank ratio
    #            --------------------------------------------------------\n")

    blanks <-  dataReduced[get(COLNAMES[["Sample_type"]]) %in% BLANK]
    Yblank <- COLNAMES[["Y"]]


#closeAllConnections()
    # prints recorded time
    startTime = Sys.time()

    dataSB <- my_fcn(
      nCORE,
      xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      y = Yorigin,
      y_trans = Y,
      yblank = Yblank,
      func = trimm_signalBlank,#(function in trimmEnds.R)
      blanks = blanks,
      noise = NOISE
    ) |> plyr::ldply(.id = NULL)

    #closeAllConnections()

    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "output", split = TRUE)
    #sink(my_log, append = TRUE, type = "message", split = TRUE)

    logr::put(paste("signal/blank: ",format(round(difftime(Sys.time(), startTime),2)), "\n"))
    Y = "Y_sb"

    processingFeature <- data.table::data.table(dplyr::full_join(dataSB, processingFeature[!IDintern %in% dataSB$IDintern], by = colnames(processingFeature)))
    processingGroup <- dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,'SignalBlank' :=any(signalBlankRatio %in% TRUE), groupIndices][,.(groupIndices, SignalBlank)]), by = c("groupIndices"))

    logr::put(paste(data.table::uniqueN(processingFeature |> dplyr::filter(signalBlankRatio %in% TRUE) %>% dplyr::select(IDintern)), " Signals were removed according to the Signal to blank ratio of ",NOISE,"."))

    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- countList[[1]]
    processingGroup <- dplyr::full_join(processingGroup, countList[[2]], by = c("groupIndices", "Feature_ID", "Batch"))

    # check length of points
    checkData <- checkLength(step, processingGroup, processingFeature, Compounds, Dilutions, Series, Signals, MIN_FEATURE)
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1

  }


  #### first outlier detection ####
  if(FOD %in% TRUE){
    logr::sep("First Outlier Detection")
    # rlang::inform("
    #               --------------------------------------------------------
    #               \tFirst Outlier Detection
    #               --------------------------------------------------------\n")
#sink()
    # prints recorded time
    startTime = Sys.time()

    dataFOD <- my_fcn(
      nCORE,
      xs = 1: data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      x = X,
      y = Y,
      func = chooseModel,#(function in FittingModel.R)
      abbr = "FOD",
      model = FOD_MODEL,
      SDRES_MIN = FOD_SDRES_MIN,
      STDRES = FOD_STDRES_MAX
    ) |> unlist(recursive = F)

    #closeAllConnections()

    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "output")
    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "message")

    logr::put(paste("FOD: ",format(round(difftime(Sys.time(), startTime),2)),"\n"))
    Y = "Y_FOD"

    dataFODModel <- tibble::tibble(
      groupIndices = as.integer(names(purrr::map(dataFOD, 1))),
      ModelName = purrr::map(dataFOD, 1) %>% unlist(use.names = F),
      Model = purrr::map(dataFOD, 2),
    )

    dataFOD = purrr::map(dataFOD,3)|> plyr::ldply(.id = NULL)

    processingFeature <- data.table::data.table(dplyr::full_join(dataFOD, processingFeature[!IDintern %in% dataFOD$IDintern], by = colnames(processingFeature)))
    processingGroup <- dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,'OutlierFOD' :=any(OutlierFOD %in% TRUE), groupIndices][,.(groupIndices, OutlierFOD)]), by = c("groupIndices"))

    logr::put(paste0("An Outlier were found for ", data.table::uniqueN(processingFeature |> dplyr::filter(OutlierFOD %in% TRUE) %>% dplyr::select(groupIndices)), " ",Series,"."))

    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- countList[[1]]
    processingGroup <- dplyr::full_join(processingGroup, countList[[2]], by = c("groupIndices", "Feature_ID", "Batch"))

    # check length of points
    checkData <- checkLength(step, processingGroup, processingFeature, Compounds, Dilutions, Series, Signals, MIN_FEATURE)
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1


  }
  #### trim ####
  if(TRIMM %in% TRUE){

    logr::sep("Trim data")
    # rlang::inform("
    #               --------------------------------------------------------
    #               \tTrim data: first Dilution should have the smallest Intensity and last point should have the biggest.
    #               --------------------------------------------------------\n")
    # #sink()
    startTime = Sys.time()
    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    dataTrim <- my_fcn(
     # cl,
      nCORE,
      xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      func = trimEnds,
      x = X,
      y = Y
    ) |>
      plyr::ldply(.id = NULL)

    #parallel::stopCluster(cl)
    #closeAllConnections()

    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "output")
    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "message")

    logr::put(paste("Trimming: ",format(round(difftime(Sys.time(), startTime),2)), "\n"))

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
    processingGroup <- dplyr::full_join(processingGroup, countList[[2]], by = c("groupIndices", "Feature_ID", "Batch"))

    logr::put(paste0("In total ", TrimFeatures," ",  Signals," in ", TrimGroups, " ",Series," were trimmed."))

    # check length of points
    checkData <- checkLength(step, processingGroup, processingFeature, Compounds, Dilutions, Series, Signals, MIN_FEATURE)
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1
    # assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))

  }

  #### second outlier detection ####

  if(SOD %in% TRUE){

    logr::sep("Second Outlier Detection")

    # rlang::inform("
    #               --------------------------------------------------------
    #               \tSecond Outlier Detection
    #               --------------------------------------------------------\n")


    processingFeature <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE, groupIndices]]
#sink()
    startTime = Sys.time()
    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    dataSOD <- my_fcn(
      #cl,
      nCORE,
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

    #closeAllConnections()

    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "output")
    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "message")
    #rm(cl)


    logr::put(paste("SOD: ",format(round(difftime(Sys.time(), startTime),2)), "\n"))

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

    logr::put(paste0("An Outlier were found for ", data.table::uniqueN(processingFeature |> dplyr::filter(OutlierSOD %in% TRUE) %>% dplyr::select(groupIndices)), " ",Series,"."))

    countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
    processingFeature <- countList[[1]]
    processingGroup <- dplyr::full_join(processingGroup, countList[[2]], by = c("groupIndices", "Feature_ID", "Batch"))

    # check length of points
    checkData <- checkLength(step, processingGroup, processingFeature, Compounds, Dilutions, Series, Signals, MIN_FEATURE)
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1

  }

  #### positive association ####
{
  logr::sep("Check if data curves have a positive slope")
    # rlang::inform("
    #               --------------------------------------------------------
    #               \tCheck if data curves have a positive slope
    #               --------------------------------------------------------\n")
  processingFeature <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE, groupIndices]]


  #sink()
    startTime = Sys.time()
    #cl <- parallel::makeCluster(getOption("cl.cores", nCORE))
    dataTrimPos <- my_fcn(
      #cl,
      nCORE,
      xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
      inputData = processingFeature,
      func = trim_pos_associated,
      x = X,
      y = Y,
      MIN_Feature = MIN_FEATURE
    ) |>
      plyr::ldply(.id = NULL)

    #closeAllConnections()

    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "output")
    #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "message")

    logr::put(paste("TrimmPos: ",format(round(difftime(Sys.time(), startTime),2)),"\n"))

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
    processingGroup <- dplyr::full_join(processingGroup, countList[[2]], by = c("groupIndices", "Feature_ID", "Batch"))

    logr::put(paste0("In total ", TrimPosFeatures," ",  Signals," in ", TrimPosGroups, " ",Series," were removed."))

    # check length of points
    checkData <- checkLength(step, processingGroup, processingFeature, Compounds, Dilutions, Series, Signals, MIN_FEATURE)
    processingFeature <- checkData[[1]]
    cutoffNew <- checkData [[2]]

    cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))

    step = step + 1
    # assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))



}




  ##### find linear Range ####
{
  logr::sep("Determining linear range")
  # rlang::inform("
  #                --------------------------------------------------------
  #                \tDetermining linear range
  #                --------------------------------------------------------\n")

  processingFeature <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE, groupIndices]]

  #dataLin <-  data.table::copy(processingFeature)[groupIndices %in% dataSODModel[aboveMinCor %in% TRUE, groupIndices]]
  #dataLin <-  data.table::copy(processingFeature)[groupIndices %in% dataSODModel[ , groupIndices]]
  #dataLin$fittingModel <- sapply(dataLin$groupIndices, function(i) dataSODModel$Model[dataSODModel$groupIndices %in% i])
#sink()

  startTime = Sys.time()
  #cl <- parallel::makeCluster(getOption("cl.cores", nCORE), type = "PSOCK")


  dataLinearRange <- my_fcn(
    #cl = myCluster,
    nCORE,
    xs = 1 : data.table::uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    func = findLinearRange,
    x = X,
    y = Y,
    real_x = COLNAMES[["X"]],
    min_feature = MIN_FEATURE,
    sd_res_factor = LR_SD_RES_FACTOR
  )

  #closeAllConnections()

  #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "output")
  #sink(paste0(OUTPUT_DIR,"/my_log.txt"), append = TRUE, type = "message")

  logr::put(paste("FindLinear Range: ",format(round(difftime(Sys.time(), startTime),2)),"\n"))

  Y = "Y_LR"

  dataLRFeature = data.table::data.table(purrr::map(dataLinearRange,1)|> plyr::ldply(.id = NULL))
  dataLRGroup = data.table::data.table(purrr::map(dataLinearRange,2)|> plyr::ldply(.id = NULL))
  dataLRGroup$aboveR2 = dataLRGroup$R2 >= R2_MIN
  dataLRGroup$LRFlag[dataLRGroup$aboveR2 %in% FALSE] = "small R2"
  dataLRGroup$enoughPointsWithinLR[dataLRGroup$aboveR2 %in% FALSE] = FALSE


  dataLRFeaturesub <- dataLRFeature[groupIndices %in% dataLRGroup[aboveR2 %in% FALSE, groupIndices] & IsLinear %in% "darkseagreen", .(IDintern, color, IsLinear, Comment)]
  dataLRFeature[dataLRFeaturesub,':=' ( color = "black",
                                        IsLinear = FALSE,
                                        Comment = paste0(Comment, "_small R2")) , on = "IDintern"]


  processingFeature <- data.table::data.table(dplyr::full_join(dataLRFeature, processingFeature[!IDintern %in% dataLRFeature$IDintern], by = colnames(processingFeature)))
  processingGroup <-  dplyr::full_join(processingGroup, dataLRGroup, by = "groupIndices")

  logr::put(paste0("For ", data.table::uniqueN(processingGroup |> dplyr::filter(aboveR2 %in% TRUE) %>% dplyr::select(groupIndices)), " ",Series,
          " a linear Range with a minimum of ", MIN_FEATURE, " Points and an R^2 higher or equal ", R2_MIN," were found."))

  countList <- countMinimumValue(processingFeature, MIN_FEATURE, step = step, y = Y)
  processingFeature <- countList[[1]]
  processingGroup <- dplyr::full_join(processingGroup, countList[[2]], by = c("groupIndices", "Feature_ID", "Batch"))



  # check length of points
  checkData <- checkLength(step, processingGroup, processingFeature, Compounds, Dilutions, Series, Signals, MIN_FEATURE)
  processingFeature <- checkData[[1]]
  cutoffNew <- checkData [[2]]

  cutoff <- dplyr::full_join(cutoff, cutoffNew, by = colnames(cutoff))
}

  ##### batch reproducibility ####
  logr::sep("check batch reproducibility")
  # rlang::inform("
  #               --------------------------------------------------------
  #               \tcheck batch reproducibility
  #               --------------------------------------------------------\n")

  {
    data.table::setDT(processingFeature)
    data.table::setDT(processingGroup)

    processingGroup <-   processingGroup[, enoughPeak_allBatches := all(get(paste0("enoughPeaks_", step)) %in% TRUE), by = Feature_ID]
    conflictBatches <- processingGroup[enoughPeak_allBatches %in% FALSE]
    processingFeature$Y_allBatches <- processingFeature$Y_LR

    if(BATCH_HARMONIZATION == TRUE){
      processingFeature[groupIndices %in% conflictBatches$groupIndices, Y_allBatches := FALSE]
      processingFeature[groupIndices %in% conflictBatches$groupIndices, color := "grey"]
      }

    processingFeature[groupIndices %in% conflictBatches$groupIndices,
                      Comment := paste(processingFeature[groupIndices %in% conflictBatches$groupIndices, Comment],"_notEnoughPeaksInAllBatches")]

    logr::put(paste0(data.table::uniqueN(conflictBatches$Feature_ID), " ", Compounds, " were excluded, because they do not show repeatable ", Series, "."))
    logr::put(paste0("Final Data set with ", data.table::uniqueN(processingGroup[enoughPeak_allBatches %in% TRUE, Feature_ID]), " ", Compounds,
            " and ", data.table::uniqueN(processingGroup[enoughPeak_allBatches %in% TRUE, Batch]), " Batch(es) -> ",
            data.table::uniqueN(processingGroup[enoughPeak_allBatches %in% TRUE, groupIndices]), " ", Series,"."))

  }

  ### Back calculation Concentration
  processingFeature  <- getLRstatus(dats =  processingFeature, datCal = processingGroup, y =  "Y")
  processingFeature <- processingFeature[!is.na(IDintern)]

  # if(TYPE %in% "targeted"){
  #   logr::put("Back calculation of concentration for the concentration series signals")
  #   Y <- ifelse(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y), "Y_trans", "Y")
  #   processingFeature$xfactor <- round(processingFeature$X/processingFeature[[COLNAMES[["X"]]]],3)
  #   processingGroup <- dplyr::full_join(processingGroup, unique(processingFeature[,c("groupIndices", "xfactor")]), by = "groupIndices")
  #   processingFeature <- getConc(dats = processingFeature, datCal = processingGroup,y = Y,INVERSE_Y =  INVERSE_Y)
  #
  #   table.backcalc.all <- processingFeature[processingFeature$Feature_ID %in% processingGroup$Feature_ID[processingGroup$enoughPeak_allBatches %in% TRUE]]
  #   table.backcalc <- table.backcalc.all[,1:9]
  #   table.backcalc$ConcentrationLR <- table.backcalc.all$ConcentrationLR
  #   table.backcalc$Status_LR <- table.backcalc.all$Status_LR
  #   table.backcalc$precision <- abs(table.backcalc$ConcentrationLR*100/table.backcalc[[COLNAMES[["X"]]]] -100)
  #
  #   t1 <- table.backcalc |> dplyr::group_by(Feature_ID, Batch) |>
  #     dplyr::summarize( lr_signals = sum(!is.na(Y)),
  #                       lr_signals_true = sum(Status_LR %in% TRUE),
  #                       'prc_all_<=20' = sum(precision <= 20)
  #     )
  #   t2 <- table.backcalc |> dplyr::group_by(Feature_ID, Batch) |>
  #     dplyr::filter(Status_LR %in% TRUE) |>
  #     dplyr::summarize( 'prc_true_<=20' = sum(precision <= 20))
  #   t <- dplyr::full_join(t1, t2, by = c("Feature_ID", "Batch") )
  #   t$prc_true_percent <- round(t$'prc_true_<=20'*100/t$lr_signals_true,2)
  #
  #
  #
  #   logr::put(paste("Concentration Back calculation", ":"))
  #   logr::put(t, n = nrow(t))
  #
  #
  #
  #
  # }





  #### biological samples ####

  if(GET_LR_STATUS %in% TRUE){ #CAL_CONC %in% TRUE |

    logr::sep("prepare biological samples")
  # rlang::inform("
  #               --------------------------------------------------------
  #               \tprepare biological samples
  #               --------------------------------------------------------\n")
    SampleFeature <- dataReduced[get(COLNAMES[["Sample_type"]]) %in% SAMPLE]

    if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)){
      SampleFeature$Y_trans <- get(TRANSFORM_Y)(SampleFeature[[column_Y]])
      SampleFeature$Y_trans[is.infinite(SampleFeature$Y_trans)] <- NA
      Y_SAMPLE <- "Y_trans"
    }

    # if(TYPE %in% "targeted" & CAL_CONC %in% TRUE){
    #   logr::put(paste("Calculate Concentration for biological Samples"))
    # SampleFeature <- getConc(dats = SampleFeature, datCal = processingGroup, y = Y_SAMPLE,INVERSE_Y =  INVERSE_Y)
    # }

    if(GET_LR_STATUS %in% TRUE){
      logr::put(" Associating biological samples to their respective linear ranges\n")
      SampleFeature  <- getLRstatus(dats = SampleFeature, datCal = processingGroup,y =  column_Y)
        #data_Sample$Status_LR <-  data.table::between(lower = data_Sample$LRStartY, x = data_Sample[, c("Area.CorrDrift")], upper = data_Sample$LREndY, NAbounds = NA)

      SampleFeature[groupIndices %in% conflictBatches$groupIndices,
                    ':=' (LRFlag = "notEnoughPeaksInAllBatches")]


    }






  }

  #### QC samples ####

  if(!is.null(QC)){

    logr::sep("check QC samples")
    # rlang::inform("
    #               --------------------------------------------------------
    #               \tcheck QC samples
    #               --------------------------------------------------------\n")


    SampleQC <- dataReduced[get(COLNAMES[["Sample_type"]]) %in% QC]

    if(!is.null(BLANK)){

      SampleBlank <- dataReduced[get(COLNAMES[["Sample_type"]]) %in% BLANK]
      SampleQC <- rbind(SampleQC, SampleBlank)
    }


    if(TRANSFORM %in% TRUE & !is.null(TRANSFORM_Y)){
      SampleQC$Y_trans <- get(TRANSFORM_Y)(SampleQC[[column_Y]])
      SampleQC$Y_trans[is.infinite(SampleQC$Y_trans)] <- NA
      Y_SAMPLE <- "Y_trans"
    }


      rsd_before <- SampleQC |>
        dplyr::group_by(Batch = get(COLNAMES[["Batch"]]), Compound = get(COLNAMES[["Feature_ID"]]), Sample.Type = get(COLNAMES[["Sample_type"]])) |>
        dplyr::summarize(.groups = "keep", rsd = sd(get(column_Y), na.rm = T)/mean(get(column_Y), na.rm = T) * 100) |>
        dplyr::group_by(Batch, Sample.Type) |>
        dplyr::summarize(.groups = "keep",median_rsd_before = median(rsd, na.rm = T))

      SampleQC  <- getLRstatus(dats = SampleQC, datCal = processingGroup,y =  column_Y)

      rsd_after <- SampleQC |>
        dplyr::filter(Status_LR %in% TRUE) |>
        dplyr::group_by(Batch = get(COLNAMES[["Batch"]]), Compound = get(COLNAMES[["Feature_ID"]]), Sample.Type = get(COLNAMES[["Sample_type"]])) |>
        dplyr::summarise(.groups = "keep",rsd = sd(get(column_Y), na.rm = T)/mean(get(column_Y), na.rm = T) * 100) |>
        dplyr::group_by(Batch, Sample.Type) |>
        dplyr::summarize(.groups = "keep",median_rsd_after = median(rsd, na.rm = TRUE))

      SampleFeature <- dplyr::full_join(SampleFeature, SampleQC, by = colnames(SampleQC))
      logr::put(paste("QC", ":"))
      logr::put(dplyr::full_join(rsd_before, rsd_after, by = c("Batch","Sample.Type")))

  }


  # assert_that(n_distinct(processList$processing$groupIndices) == rawCompounds)
  #### output files ####
  logr::sep("prepare output files")
  # rlang::inform("
  #               --------------------------------------------------------
  #               \tprepare output files
  #               --------------------------------------------------------\n")

  cutoff <- dplyr::full_join(cutoff,
                             processingFeature,
                             by = colnames(cutoff))


  processingFeature <- dplyr::full_join(processingFeature, cutoff, by = colnames(cutoff)) |>
    tidyr::drop_na(IDintern)


  assertthat::are_equal(nrow(processingFeature), nrow(processingFeatureCal))

  #output1.1 <- data.table::copy(processingFeature)
  #output2.1 <- data.table::copy(processingGroup)
  #output3.1 <- data.table::copy(SampleFeature)


  #output:
  data.table::setnames(processingFeature, c("Feature_ID","Sample_ID","Class","Sample.Type", "Batch", "Y", "Injection_order"), c(COLNAMES[["Feature_ID"]],COLNAMES[["Sample_ID"]],COLNAMES[["Class"]],COLNAMES[["Sample_type"]], COLNAMES[["Batch"]], COLNAMES[["Y"]], COLNAMES[["Injection_order"]] ))
  data.table::setnames(processingGroup, c("Feature_ID", "Batch"), c(COLNAMES[["Feature_ID"]], COLNAMES[["Batch"]]))


  #1)full table dilution/concentration curves - Signal based
  output1 <- data.table::copy(processingFeature)
  #output1.1 <- processingFeature

  #2) full table dilution/concentration curves - Feature based
  output2 <- data.table::copy(processingGroup)
  #output2.1 <- processingGroup


  # #3) filtered table dilution/concentration curves - Signal based (high quality)
  # output3 <- processingFeature[groupIndices %in% processingGroup[
  #   get(paste0("enoughPeaks_",step)) %in% TRUE, groupIndices]]
  # output3 <- output3 |> dplyr::group_by(groupIndices) |> dplyr::filter(!all(is.na(Y_allBatches)))
  #

  #4) full table biological Samples - Signal based

  if(!is.null(column_sampleClass))dataOrigin[[column_sampleClass]] <- as.character(dataOrigin[[column_sampleClass]])

  table_Feature_all <- dplyr::full_join(data.table::copy(SampleFeature),data.table::copy(processingFeature), by = colnames(SampleFeature))
  output3 <- dplyr::full_join(data.table::copy(table_Feature_all), dataOrigin, by = intersect(colnames(dataOrigin), colnames(table_Feature_all)))

  output4 <- output3[,c(colnames(dataOrigin), "groupIndices","Status_LR", "LRFlag"), with = F]
  #output3.1 <- SampleFeature

  # #5) filtered table biological Samples - Signal based (high quality)
  # output5 <- SampleFeature |> dplyr::filter(get(COLNAMES[["Feature_ID"]]) %in% unlist(unique(output3[COLNAMES[["Feature_ID"]]])))

  #6) summary table
  # output4 <- data.table::copy(SampleFeature) |>
  #   #subset(get(COLNAMES[["Sample_type"]]) %in% QC) |>
  #   dplyr::group_by(Feature_ID = get(COLNAMES[["Feature_ID"]]),
  #                   Batch = get(COLNAMES[["Batch"]]),
  #                   Type = get(COLNAMES[["Sample_type"]]),
  #                   Class = get(COLNAMES[["Class"]]),
  #                   LRFlag) |>
  #   dplyr::summarize(.groups = "drop",
  #     'Signals' = length(Status_LR),
  #     'LR_Status' = vctrs::vec_count(Status_LR)$key,
  #     'LR_Status_n' = vctrs::vec_count(Status_LR)$count,
  #     'LR_Status_n[%]' = round(vctrs::vec_count(Status_LR)$count/length(Status_LR)*100,2),
  #     rsd_all = round(sd(get(Y_SAMPLE), na.rm = T)/mean(get(Y_SAMPLE), na.rm = T) * 100,2),
  #     rsd_LR_TRUE = round(sd(get(Y_SAMPLE)[Status_LR %in% TRUE], na.rm = T)/mean(get(Y_SAMPLE)[Status_LR %in% TRUE], na.rm = T) * 100,2),
  #    ) #|>
    #
    # tidyr::pivot_wider(names_from = c(Type),
    #                    values_from = c('LR_TRUE', 'LR_TRUE[%]', rsd_all, rsd_LR_TRUE)
    # ) |>
    #
    # dplyr::select(tidyr::contains(unique(SampleFeature$Sample.Type)) )


  # htmloutput6 <-  DT::datatable(output6) %>%
  #   DT::formatStyle(dplyr::select(output6,tidyr::contains("%"))|> colnames(),
  #               backgroundColor = DT::styleInterval(c(20,80), c('red','yellow', 'green'))
  #   )
  # htmlwidgets::saveWidget(htmloutput6,
  #                         file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "Compounds_summary.html" , sep = "_")),
  #                         selfcontained = TRUE)
  #
  # output6.1 <- SampleFeature |>
  #   subset(get(COLNAMES[["Sample_type"]]) %in% SAMPLE ) |>
  #   dplyr::group_by(Feature_ID = get(COLNAMES[["Feature_ID"]]), Batch = get(COLNAMES[["Batch"]])) |>
  #   dplyr::select(Sample_ID = all_of(SAMPLE_ID), Y = all_of(Y_SAMPLE)) |>
  #   tidyr::pivot_wider(names_from = Sample_ID,
  #                      values_from = Y
  #   ) |>
  #   dplyr::full_join(output6 |> dplyr::select(Feature_ID, Batch, 'LR_TRUE[%]_Sample'), . , by = c("Feature_ID", "Batch")) |>
  #   dplyr::select("Feature_ID", "Batch", "LR_TRUE[%]_Sample", everything())
  #
  # htmloutput6.1 <-  DT::datatable(output6.1, filter = 'top',
  #                             extensions = 'Buttons', options = list(
  #                               pageLength = nrow(output6.1),
  #                               dom = 'Bfrtip',
  #                               buttons = c('copy', 'csv', 'excel', 'print')
  #                             )) %>%
  #   DT::formatStyle("LR_TRUE[%]_Sample",
  #               backgroundColor = DT::styleInterval(c(20,80), c('red','yellow', 'green'))
  #   )
  # htmlwidgets::saveWidget(htmloutput6.1, file.path( REPORT_OUTPUT_DIR, paste(Sys.Date(), PREFIX, "Compounds_biol_samples.html" , sep = "_")), selfcontained = FALSE)





  data.table::setnames(skip_absent = T, processingGroup, c("Feature_ID","Sample_ID", "Batch", "Y", "X", "Injection_order"),
                       c(COLNAMES[["Feature_ID"]],COLNAMES[["Sample_ID"]], COLNAMES[["Batch"]], COLNAMES[["Y"]], COLNAMES[["X"]], COLNAMES[["Injection_order"]] ))

  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  if(any(colnames(processingFeature) %in% "Y_trans")){
    colnames(output1)[which(colnames(output1) %in% "Y_trans")] <- paste0("Y_transformed(",TRANSFORM_Y, ")")
    colnames(output2)[which(colnames(output2) %in% "Y_trans")] <- paste0("Y_transformed(",TRANSFORM_Y, ")")
    colnames(output3)[which(colnames(output3) %in% "Y_trans")] <- paste0("Y_transformed(",TRANSFORM_Y, ")")
    #colnames(output4)[which(colnames(output4) %in% "Y_trans")] <- paste0("Y_transformed(",TRANSFORM_Y, ")")
    colnames(SampleQC)[which(colnames(SampleQC) %in% "Y_trans")] <- paste0("Y_transformed(",TRANSFORM_Y, ")")



  }

  if(any(colnames(processingFeature) %in% "X_trans")){
    colnames(output1)[which(colnames(output1) %in% "X_trans")] <- paste0("X_transformed(",TRANSFORM_X, ")")
    colnames(output2)[which(colnames(output2) %in% "X_trans")] <- paste0("X_transformed(",TRANSFORM_X, ")")
    colnames(output3)[which(colnames(output3) %in% "X_trans")] <- paste0("X_transformed(",TRANSFORM_X, ")")
    #colnames(output4)[which(colnames(output4) %in% "X_trans")] <- paste0("X_transformed(",TRANSFORM_X, ")")
    colnames(SampleQC)[which(colnames(SampleQC) %in% "X_trans")] <- paste0("X_transformed(",TRANSFORM_X, ")")

  }


  if(GET_OUTPUT %in% TRUE){

    if(any(c("DilutionCurves", "all") %in% which_output)){
      write.csv(output1, file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", Series,"_signalBased.csv")))
      write.csv(output2, file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", Series,"_featureBased.csv")))
      # writexl::write_xlsx(x = list(Signals = output1, Features = output2),
      #                     path = file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", Series,".xlsx")))
    }

    if(any(c("BiologicalSamples", "all") %in% which_output)){
     # write.csv(output3, file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", "BiologicalSamples_signalBased.csv")))
      write.csv(output4, file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", "result.csv")))

      # writexl::write_xlsx(x = list(Signals = output3, summary = output4),
      #                     path = file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", "BiologicalSamples.xlsx")))
    }

  }

  ## Plots

  #7) bar plot summary per dilution/concentration

  printPlot <- ifelse(any(c("Plots", "all") %in% which_output), TRUE, FALSE)
  if(is.null(TRANSFORM_Y)){
    Yraw <- COLNAMES[["Y"]]
  } else {
    Yraw <- paste0("Y_transformed(",TRANSFORM_Y, ")")
  }

  if(is.null(TRANSFORM_X)){
    Xraw <- COLNAMES[["X"]]
  } else {
    Xraw <- paste0("X_transformed(",TRANSFORM_X, ")")
  }



  summary_barplot <- plot_Barplot_Summary(printPDF = printPlot,
                                          inputData_Series = output1,
                                          COLNAMES = COLNAMES,
                                          X = Xraw, Y = Yraw,
                                          output_dir = IMG_OUTPUT_DIR,outputfileName = paste0(PREFIX,"_Summary_Calibration_Barplot"))
  logr::put("summary_barplot was created")

  #8) scatter plot
  FDS_scatterplot <- plot_FDS(printPDF = printPlot,
                              inputData_Series = output1,
                              inputData_BioSamples = output3 |> dplyr::filter(get(COLNAMES[["Sample_type"]]) %in% SAMPLE),
                              inputData_QC = SampleQC,
                              COLNAMES = COLNAMES, Xcol = Xraw, Ycol = Yraw, TRANSFORM_Y = TRANSFORM_Y, inverse_y = INVERSE_Y,
                              Series = Series, output_dir = IMG_OUTPUT_DIR,outputfileName = paste0(PREFIX,"_CalibrationPlot")
  )

  logr::put("Scatterplot was created")
  #9) bar plot summary for all samples


  summary_barplot_all <- plot_Barplot_Summary_Sample(printPDF = printPlot,
                                                     inputData_Samples = output3,
                                                     COLNAMES = COLNAMES,
                                                     X = Xraw, Y = Yraw,
                                                     output_dir = IMG_OUTPUT_DIR,
                                                     group = "Batch",
                                                     ordered = column_injectionOrder,
                                                     outputfileName = paste0(PREFIX,"_Summary_Barplot_All"))

  logr::put("summary_barplot_all was created")

  #10) bar plot summary for biological samples

  summary_barplot_sample <- plot_Barplot_Summary_Sample(printPDF = printPlot,
                                                        inputData_Samples = output3 |> dplyr::filter(get(COLNAMES[["Sample_type"]]) %in% SAMPLE),
                                                        COLNAMES = COLNAMES,
                                                        X = Xraw, Y = Yraw,
                                                        output_dir = IMG_OUTPUT_DIR,
                                                        group = "Batch",
                                                        Plotfill = "Class",
                                                        ordered = "Class",
                                                        outputfileName = paste0(PREFIX,"_Summary_Barplot_Samples"))

  logr::put("summary_barplot_sample was created")


  #11) barplot summary for QC samples

  summary_barplot_QC <- plot_Barplot_Summary_Sample(printPDF = printPlot,
                                                    inputData_Samples = SampleQC,
                                                    COLNAMES = COLNAMES,
                                                    X = Xraw, Y = Yraw,
                                                    output_dir = IMG_OUTPUT_DIR,
                                                    group = "Batch",
                                                    Plotfill = "Sample.Type",
                                                    ordered = "Sample.Type",
                                                    outputfileName = paste0( PREFIX,"_Summary_Barplot_QC"))

  logr::put("summary_barplot_sample was created")





  processList <- list(
    "All_DilutionCurves_Signals" = output1,
    "All_DilutionCurves_Features" = output2,
    #"All_Samples_Signals" = output3,
    "result" = output4,
    "dataModel_FOD" = dataFODModel,
    "dataModel_SOD" = dataSODModel#,
    # "Parameters" = list(
    #   data = DAT,
    #   COLNAMES = COLNAMES,
    #   nCORE = nCORE,
    #   LOG_TRANSFORM = LOG_TRANSFORM,
    #   R2min = R2min,
    #   res = res,
    #   #calcFor = calcFor,
    #   #SRES = SRES,
    #   #R2FOD = R2FOD,
    #   #R2SOD = R2SOD,
    #   MIN_FEATURE = MIN_FEATURE,
    #   X = X,
    #   Y = Y
    # )
  )

  save(processList, file = file.path( REPORT_OUTPUT_DIR, paste0(Sys.Date(),"_", PREFIX,"_", Series,"_image.RData")))

  #sink()
  logr::log_close()

  return(processList)

}

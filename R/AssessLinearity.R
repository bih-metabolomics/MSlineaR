#' Title
#'
#' @param COLNAMES
#' @param DAT data set with minimum following columns must be provided: ID, REPLICATE,Concentration or Dilution and INtensity or Area
#' @param NCORE number of cores to use for parallelization. Default value is 1
#' @param MIN_FEATURE minimum number of features within a dilution or concentration series to be considered as linear range
#' @param LOG_TRANSFORM logical Value. should the data be log transformed? Default is TRUE.
#' @return
#' @export
#'
#' @examples
AssessLinearity <- function(COLNAMES = c(ID =  "featNames",
                                         #Batch = "Batch",
                                         REPLICATE = NULL,
                                         X = "Concentration",
                                         Y =   "area.raw"),
                            DAT,
                            nCORE = 1,
                            MIN_FEATURE = 5,
                            LOG_TRANSFORM = TRUE,
                            R2min = 0.95,
                            ...) {
  devtools::load_all()

  source(file = "R/prep.R")
  source(file = "R/countMinimumValues.R")
  source(file = "R/findLinearRange.R")
  source(file = "R/FittingModel.R")
  # source(file = "R/getSummary.R")
  # source(file = "R/outlierDetection.R")
  # source(file = "R/plotSignals.R")
  source(file = "R/trimEnds.R")

  dataOrigin <-checkData(data.table::data.table(DAT), MIN_FEATURE, LOG_TRANSFORM, nCORE)

  if(LOG_TRANSFORM %in% TRUE){
    X = "XLog"
    Y = "YLog"
  } else{
    X = "X"
    Y = "Y"
  }

  processingFeature <- data.table::copy(dataOrigin)

  #### normalizing, centralizing, log transforming ####
   message("preparing data\n--------------------------------------------------------\n")

  dataPrep <- prepareData(processingFeature)

  processingFeature <- data.table::copy(dataPrep)

  countList <- countMinimumValue(processingFeature, MIN_FEATURE)
  processingFeature <- data.table::data.table(countList[[1]])
  processingGroup <- data.table::data.table(countList[[2]])


   # save key numbers
  nCompounds <- data.table::uniqueN(processingGroup, by = c("ID"))
  nReplicates <- data.table::uniqueN(processingFeature, by = c("REPLICATE"))
  nDilutions <- data.table::uniqueN(processingFeature[, DilutionPoint := 1:.N, by = c("ID", "REPLICATE")], by = c("DilutionPoint"))
  nSeries <- data.table::uniqueN(processingGroup, by = c("groupIndices"))
  nPeaks <- data.table::uniqueN(processingFeature, by = c("IDintern"))

  message(
    "Data set with ", nCompounds, " molecular Features, ",
    nReplicates, " REPLICATE(s) and ",
    nDilutions, " X / Dilutions -> ",
    format.default(nSeries, big.mark = ",", scientific = F), " Feature dilution series.",
    "\n--------------------------------------------------------\n"
  )

  nCompoundsNew <- data.table::uniqueN(processingGroup[enoughPeaks %in% TRUE], by = c("ID"))
  nReplicatesNew <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("REPLICATE"))
  nSeriesNew <- data.table::uniqueN(processingGroup[enoughPeaks %in% TRUE], by = c("groupIndices"))

  message("Removed ", nCompounds - nCompoundsNew, " molecular Features with less than ", MIN_FEATURE, " points.
          New Data set with ", nCompoundsNew, " molecular Features and ", nReplicatesNew, " Replicates -> ", nSeriesNew, " Feature dilution series \n--------------------------------------------------------\n")

  stopifnot(exprs = {
    "all Compounds were removed" = nCompoundsNew - data.table::uniqueN(processingGroup[enoughPeaks %in% FALSE], by = c("ID")) > 0
    "all Dilution/Concentration-Series were removed" = nSeriesNew - data.table::uniqueN(processingGroup[enoughPeaks %in% FALSE], by = c("groupIndices")) > 0
  })


  # assert_that(nrow(processList$processing) == rawCompounds*n_dilution*nReplications)

  #### outlier ####

  message("First Outlier Detection\n")

  # pblapply(1:uniqueN(processing, by = "groupIndices"), function(i) processing[groupIndices %in% i, outlierDetection(.SD)]) |> ldply()

  if (nCORE > 1) future::plan(multisession, workers = nCORE)
  dataFOD <- my_fcn(
    xs = 1:data.table::uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    x = X,
    y = Y,
    func = chooseModel,
    abbr = "FOD"
  ) |> unlist(recursive = F)

  if (nCORE > 1) future::plan(sequential)


  dataFODModel <- tibble::tibble(
    groupIndices = as.integer(names(map(dataFOD, 1))),
    ModelName = map(dataFOD, 1) %>% unlist(use.names = F),
    Model = map(dataFOD, 2),
    #RMSE = map(dataFOD,3) %>% unlist(use.names = T),
    R2 = map(dataFOD, 3) %>% unlist(use.names = F)
    )

  dataFOD = map(dataFOD,5)|> plyr::ldply(.id = NULL)

  processingFeature <- data.table::data.table(dplyr::full_join(dataFOD, processingFeature[!IDintern %in% dataFOD$IDintern]))
  processingGroup <- dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,OutlierFOD :=any(OutlierFOD %in% TRUE), groupIndices][,.(groupIndices, OutlierFOD)]))

  message("An Outlier were found for ", data.table::uniqueN(processingFeature |> dplyr::filter(OutlierFOD %in% TRUE) %>% dplyr::select(groupIndices)), " Compounds.\n")

  countList <- countMinimumValue(processingFeature, MIN_FEATURE)
  processingFeature <- countList[[1]]
  processingGroup <- dplyr::full_join(subset(processingGroup[groupIndices %in% countList[[2]]$groupIndices], select= -c(N, enoughPeaks)), countList[[2]])

  #### trim ####
  message("Trim data: first Dilution should have the smallest Intensity and last point should have the biggest.")

  if (nCORE > 1) future::plan(multisession, workers = nCORE)
  dataTrim <- my_fcn(
    xs = 1: data.table::uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    func = trimEnds,
    x = X,
    y = Y
  ) |>
    plyr::ldply(.id = NULL)

  if (nCORE > 1) future::plan(sequential)

  # setDT(dataTrim)[Y < levelnoise, color := "grey"]
  # setDT(dataTrim)[Y < levelnoise, Comment := "below noise level"]


  processingFeature <- data.table::data.table(dplyr::full_join(dataTrim, processingFeature[!IDintern %in% dataTrim$IDintern]))
  processingGroup <- dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,trim :=any(trim %in% TRUE),groupIndices][,.(groupIndices, trim)]))


  # check length of points
  countList <- countMinimumValue(processingFeature, MIN_FEATURE )
  processingFeature <- dplyr::full_join(countList[[1]], processingFeature[!IDintern %in% countList[[1]]$IDintern])
  processingGroup <- dplyr::full_join(subset(processingGroup[groupIndices %in% countList[[2]]$groupIndices], select= -c(N, enoughPeaks)), countList[[2]])


  # assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))



  TrimFeatures <- data.table::uniqueN(dataTrim |> dplyr::filter(!color %in% "black") %>% dplyr::select(groupIndices))
  message(TrimFeatures, " individual Features were removed after trimming the data")


  discardCompound <- length(processingGroup[enoughPeaks %in% FALSE, groupIndices])
  remainingCompound <- length(processingGroup[enoughPeaks %in% TRUE, groupIndices])
  message(discardCompound, " Compounds had less than ", MIN_FEATURE, " Peaks and were discarded.\n--------------------------------------------------------\n")
  message("Data set with ", remainingCompound, " Compounds.\n--------------------------------------------------------\n")



  stopifnot(exprs = {
    "all Compounds were removed" = remainingCompound > 0
  })




  #### Fitting ####

  # 1.) choose Model according to correlation threshold

  message("Fitting linear, logistic and quadratic regression\n")

  if (nCORE > 1) future::plan(multisession, workers = nCORE)
  dataSOD <- my_fcn(
    xs = 1:data.table::uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    x = X,
    y = Y,
    func = chooseModel,
    abbr = "SOD",
    model = "linear",
    R2min = 0.95
  ) |> unlist(recursive = F)

  if (nCORE > 1) future::plan(sequential)


  dataSODModel <- tibble::tibble(
    groupIndices = as.integer(names(map(dataSOD, 1))),
    ModelName = map(dataSOD, 1) %>% unlist(use.names = F),
    Model = map(dataSOD, 2),
    R2 = map(dataSOD, 3) %>% unlist(use.names = F),
    aboveMinCor = map(dataSOD, 4) %>% unlist(use.names = F)

  )
  data.table::setDT(dataSODModel)
  dataSOD = map(dataSOD,5)|> plyr::ldply(.id = NULL)

  processingFeature <- data.table::data.table(dplyr::full_join(dataSOD, processingFeature[!IDintern %in% dataSOD$IDintern]))
  processingGroup <-  dplyr::full_join(processingGroup, unique(data.table::copy(processingFeature)[,OutlierSOD :=any(OutlierSOD %in% TRUE),groupIndices][,.(groupIndices, OutlierSOD)]))
  processingGroup <-  dplyr::full_join( processingGroup, subset(dataSODModel, select = -c(Model)), by = "groupIndices")


  # check length of points
  countList <- countMinimumValue(processingFeature, MIN_FEATURE )
  processingFeature <- dplyr::full_join(countList[[1]], processingFeature[!IDintern %in% countList[[1]]$IDintern])
  processingGroup <- dplyr::full_join(subset(processingGroup[groupIndices %in% countList[[2]]$groupIndices], select= -c(N, enoughPeaks)), countList[[2]])


  discardCompoundFitting <- data.table::uniqueN(processingGroup[aboveMinCor %in% FALSE, groupIndices])
  remainingCompoundFitting <- data.table::uniqueN(processingGroup[aboveMinCor %in% TRUE, groupIndices])

  message(discardCompoundFitting, " Compounds have a low R^2 (less than ", R2min, ") and will be excluded.\n")


    #### Combining ####
    # discardCompoundFitting <- dplyr::n_distinct(processing |> filter(aboveCorFit == FALSE) %>% dplyr::select(groupIndices))
    # savedCompounds <- n_distinct(processing |> filter(aboveCorFit == TRUE, enoughPeaks %in% TRUE) %>% dplyr::select(groupIndices))
    # message(savedCompounds, " of ", discardCompoundFitting, " Features have a R^2 above ", R2SOD, " after second outlier detection\n--------------------------------------------------------\n")

  #   # combine data for both outlier detections
  #   dataModelCombined <- copy(dataModel)
  #   if (n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices]) > 0) {
  #     dataModelCombined <- dataModelCombined[
  #       groupIndices %in% dataModelSOD[aboveMinCor %in% TRUE, groupIndices],
  #       ":="(Model = dataModelSOD[aboveMinCor %in% TRUE, Model],
  #         correlation = dataModelSOD[aboveMinCor %in% TRUE, correlation],
  #         aboveMinCor = dataModelSOD[aboveMinCor %in% TRUE, aboveMinCor],
  #         fittingModel = dataModelSOD[aboveMinCor %in% TRUE, fittingModel])
  #     ]
  #   }
  # } else {
  #   processing$outlierSOD <- NA
  #   dataSOD <- NULL
  #   dataConsSOD <- NULL
  #   dataModelSOD <- NULL
  #   dataModelCombined <- copy(dataModel)
  # }
  # assert_that(n_distinct(processList$processing |> filter(aboveMinCor %in% TRUE) |> dplyr::select(groupIndices)) == remainingCompoundFitting + savedCompounds )

  message("Data set with ", data.table::uniqueN(processingGroup[aboveMinCor %in% TRUE, groupIndices]), " Compounds with a R^2 above ", R2min, "\n--------------------------------------------------------\n")

  ## find linear Range
  message("Determining linear range\n--------------------------------------------------------\n")

  dataLin <-  data.table::copy(processingFeature)[groupIndices %in% dataSODModel[aboveMinCor %in% TRUE, groupIndices]]
  dataLin$fittingModel <- sapply(dataLin$groupIndices, function(i) dataSODModel$Model[dataSODModel$groupIndices %in% i])

  if (nCORE > 1) future::plan(multisession, workers = nCORE)
  dataLinearRange <- my_fcn(
    xs = 1:data.table::uniqueN(dataLin$groupIndices),
    inputData = dataLin,
    func = findLinearRange,
    x = X,
    y = Y,
    modelObject = "fittingModel",
    MIN_FEATURE = MIN_FEATURE
  ) #|>
    #unlist(recursive = F)
    #plyr::ldply(.id = NULL)

  if (nCORE > 1) future::plan(sequential)


  dataLR = map(dataLinearRange,1)|> plyr::ldply(.id = NULL)
  processingFeature <- data.table::data.table(dplyr::full_join(dataLR, processingFeature[!IDintern %in% dataLR$IDintern]))
  processingGroup <-  dplyr::full_join(processingGroup, map(dataLinearRange,2)|> plyr::ldply(.id = NULL), by = "groupIndices")





  # assert_that(n_distinct(processList$processing$groupIndices) == rawCompounds)

  message("For ", data.table::uniqueN(processingGroup[enoughPointsWithinLR %in% TRUE,groupIndices]), " Features a linear Range with a minimum of ", MIN_FEATURE, " Points were found.\n--------------------------------------------------------\n")

  #
  # <!-- plotSignals(DAT = processList$Preprocessed$dataLinearRange %>% filter(ID %in% metabolitesRandomSample), x = "DilutionPoint", y = "IntensityNorm") -->
  #
  message("summarize all informations in one list\n--------------------------------------------------------\n")

 # processing <- full_join(processing, dataPrep[!groupIndices %in% dropNA$groupIndices, -c("pch", "color", "N", "Comment", "enoughPeaks")],
#    by = c("groupIndices", "IDintern", "Intensity", "X", "IntensityNorm", "DilutionPoint", "IntensityLog", "ConcentrationLog")
#  )


 # processing <- full_join(x = processingFeature, y = dataPrep[groupIndices %in% dropNA$groupIndices, -c("N")], by = names("y"))


  data.table::setorder(processingFeature, DilutionPoint)
  dataSummary <- getSummaryList(processing)

  processing <- dplyr::full_join(
    x = dataOrigin, y = processingFeature,
    by = setNames(
      c("ID", "IDintern", "REPLICATE", "X", "Y"),
      c(COLNAMES[["ID"]], "IDintern", COLNAMES[["REPLICATE"]], COLNAMES[["X"]], COLNAMES[["Y"]])
    ),
    suffix = c(".raw", ".processed")
  )


  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  processList <- list(
    "dataOrigin_1" = dataOrigin,
    "dataPrep_2" = dataPrep,
    "dataFOD_3" <- dataFOD,
    "dataTrim_4" = dataTrim,
    "dataSOD_5" = dataSOD,
    "dataModel_6" = dataSODModel,
    "dataLinearRange_7" = dataLinearRange,
    "dataLR_8" = dataLR,
    "summaryFFDS" = processingFeature,
    "summaryFDS" = processingGroup,
    "Parameters" = list(
      data = DAT,
      COLNAMES = COLNAMES,
      nCORE = nCORE,
      #calcFor = calcFor,
      #SRES = SRES,
      #R2FOD = R2FOD,
      #R2SOD = R2SOD,
      MIN_FEATURE = MIN_FEATURE,
      X = X,
      Y = Y
    )
  )


  #if (!empty(dataFOD)) processList$"dataFOD_3" <- dataFOD
  #if (!is.null(dataSOD) | !empty(dataSOD)) processList$"dataSOD_7" <- dataSOD
  #if (!is.null(dataConsSOD) | !empty(dataConsSOD)) processList$"dataConsSOD_8" <- dataConsSOD
  #if (!is.null(dataModelSOD) | !empty(dataModelSOD)) processList$"dataModelSOD_9" <- dataModelSOD


  return(processList)
}

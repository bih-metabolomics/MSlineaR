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
                            ...) {
  devtools::load_all()

  # source(file = "R/countMinimumValues.R")
  # source(file = "R/findLinearRange.R")
  # source(file = "R/FittingModel.R")
  # source(file = "R/getSummary.R")
  # source(file = "R/outlierDetection.R")
  # source(file = "R/plotSignals.R")
  # source(file = "R/trimEnds.R")

  dataOrigin <-checkData(data.table::data.table(DAT), MIN_FEATURE, LOG_TRANSFORM, nCORE)

  if(LOG_TRANSFORM %in% TRUE){
    X = "XLog"
    Y = "YLog"
  } else{
    X = "X"
    Y = "Y"
  }

  processingFeature <- copy(dataOrigin)

  ## normalizing, centralizing, log transforming
   message("preparing data\n--------------------------------------------------------\n")

  processingFeature <- prepareData(processingFeature)
  processingGroup <- countMinimumValue(processingFeature)

  processingFeature$color[processingFeature$groupIndices %in% processingGroup$groupIndices[processingGroup$enoughPeaks %in% FALSE]] <- "grey"
  processingFeature$Comment[processingFeature$color %in% "black"] <- unlist(apply(cbind(processingFeature$Comment, "EnoughPeaks"), 1, function(x) paste(x[!is.na(x)], collapse = " _ ")))
  processingFeature$Comment[processingFeature$color %in% "grey"] <- unlist(apply(cbind(processingFeature$Comment, "notEnoughPeaks"), 1, function(x) paste(x[!is.na(x)], collapse = " _ ")))


   # save key numbers
  nCompounds <- uniqueN(processingGroup, by = c("ID"))
  nReplicates <- uniqueN(processingFeature, by = c("REPLICATE"))
  nDilutions <- uniqueN(processingFeature[, DilutionPoint := 1:.N, by = c("ID", "REPLICATE")], by = c("DilutionPoint"))
  nSeries <- uniqueN(processingGroup, by = c("groupIndices"))
  nPeaks <- uniqueN(processingFeature, by = c("IDintern"))

  message(
    "Data set with ", nCompounds, " molecular Features, ",
    nReplicates, " REPLICATE(s) and ",
    nDilutions, " X / Dilutions -> ",
    format.default(nSeries, big.mark = ",", scientific = F), " Feature dilution series.",
    "\n--------------------------------------------------------\n"
  )

  nCompoundsNew <- uniqueN(processingGroup[enoughPeaks %in% TRUE], by = c("ID"))
  nReplicatesNew <- uniqueN(processingFeature[color %in% "black"], by = c("REPLICATE"))
  nSeriesNew <- uniqueN(processingGroup[enoughPeaks %in% TRUE], by = c("groupIndices"))

  message("Removed ", nCompounds - nCompoundsNew, " molecular Features with less than ", MIN_FEATURE, " points.
          New Data set with ", nCompoundsNew, " molecular Features and ", nReplicatesNew, " Replicates -> ", nSeriesNew, " Feature dilution series \n--------------------------------------------------------\n")

  stopifnot(exprs = {
    "all Compounds were removed" = nCompoundsNew - uniqueN(processingGroup[enoughPeaks %in% FALSE], by = c("ID")) > 0
    "all Dilution/Concentration-Series were removed" = nSeriesNew - uniqueN(processingGroup[enoughPeaks %in% FALSE], by = c("groupIndices")) > 0
  })


  # assert_that(nrow(processList$processing) == rawCompounds*n_dilution*nReplications)

  ## outlier

  message("First Outlier Detection\n")

  # pblapply(1:uniqueN(processing, by = "groupIndices"), function(i) processing[groupIndices %in% i, outlierDetection(.SD)]) |> ldply()

  if (nCORE > 1) plan(multisession, workers = nCORE)
  dataFOD <- my_fcn(
    xs = 1:data.table::uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    x = X,
    y = Y,
    func = outlierDetection,
    numboutlier = 1,
    SRES = 2
  ) |>
    plyr::ldply(.id = NULL)

  plan(sequential)
  setDT(dataFOD)[outlier %in% TRUE, Comment := str_replace(Comment, "outlier$", "outlierFOD")]
  # assert_that(n_distinct(dataFOD$groupIndices) == n_distinct(processing$groupIndices))
  processingFeature[, Comment := as.character(Comment)][IDintern %in% dataFOD$IDintern, ":="(Comment = dataFOD$Comment,
    color = dataFOD$color,
    pch = dataFOD$pch,
    outlierFOD = dataFOD$outlier), on = "IDintern"]


  message("An Outlier were found for ", uniqueN(processingFeature |> dplyr::filter(outlierFOD %in% TRUE) %>% dplyr::select(groupIndices)), " Compounds.\n")

  ## trim
  message("Trim data: first Dilution should have the smallest Intensity and last point should have the biggest.")

  if (nCORE > 1) plan(multisession, workers = nCORE)
  dataTrim <- my_fcn(
    xs = 1:uniqueN(processingFeature$groupIndices),
    inputData = processingFeature,
    func = trimEnds,
    x = X,
    y = Y
  ) |>
    plyr::ldply(.id = NULL)

  setDT(dataTrim)[Intensity < levelnoise, color := "grey"]
  setDT(dataTrim)[Intensity < levelnoise, Comment := "below noise level"]


  plan(sequential)

  # assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))

  processing[, Comment := as.character(Comment)][IDintern %in% dataTrim$IDintern, ":="(Comment = dataTrim$Comment,
    color = dataTrim$color,
    pch = dataTrim$pch,
    outlierFOD = dataTrim$outlierFOD,
    trim = dataTrim$trim), on = "IDintern"]

  # check length of consecutive points
  message("check length of consecutive points")
  TrimFeatures <- n_distinct(dataTrim |> dplyr::filter(!color %in% "black") %>% dplyr::select(groupIndices))
  message(TrimFeatures, "individual Features were removed after trimming the data")

  dataCons <- processing[color == "black", .N, by = .(groupIndices)]
  dataCons[, enoughPeaks := N >= MIN_FEATURE][enoughPeaks %in% F, Comment := "notEnoughPeaks"]

  # assert_that(n_distinct(dataCons$groupIndices) == n_distinct(processing$groupIndices))

  processing[
    groupIndices %in% dataCons$groupIndices,
    enoughPeaks := rep(dataCons$enoughPeaks, each = nDilutions)
  ]

  processing[
    groupIndices %in% dataCons$groupIndices,
    CommentTrim := rep(dataCons$Comment, each = nDilutions)
  ]

  processing[Comment %in% c(NA, NULL, "", " "), Comment := NA]
  processing$Comment <- unlist(apply(cbind(processing$Comment, processing$CommentTrim), 1, function(x) paste(x[!is.na(x)], collapse = " _ ")))




  # processing[ , Comment := as.character(Comment)][groupIndices %in% dataCons$groupIndices & (!Comment %in% c(NA, NULL, "", " ")), Comment := paste(Comment,dataCons$Comment, sep = "_")]
  # processing[ , Comment := as.character(Comment)][groupIndices %in% dataCons$groupIndices & enoughPeaks %in% FALSE & (Comment %in% c(NA, NULL, "", " ")), Comment := "notEnoughPeaks"]




  discardCompound <- n_distinct(processing %>% filter(enoughPeaks == FALSE, !color %in% "black") %>% dplyr::select(groupIndices))
  remainingCompound <- n_distinct(processing %>% filter(enoughPeaks == TRUE, color %in% "black") %>% dplyr::select(groupIndices))
  message(discardCompound, " Compounds had less than ", MIN_FEATURE, " Peaks and were discarded.\n--------------------------------------------------------\n")
  message("Data set with ", remainingCompound, " Compounds.\n--------------------------------------------------------\n")

  ## Fitting

  # 1.) choose Model according to correlation threshold

  message("Fitting linear, logistic and quadratic regression\n")

  if (nCORE > 1) plan(multisession, workers = nCORE)
  dataModel <- my_fcn(
    xs = 1:n_distinct(processing |> filter(enoughPeaks %in% TRUE, color %in% "black") |> dplyr::select(groupIndices)),
    inputData = processing[enoughPeaks %in% TRUE & color %in% "black"],
    func = chooseModel,
    x = X,
    y = Intensity
  ) %>% unlist(recursive = F)

  plan(sequential)

  # gc()

  dataModel <- tibble(
    groupIndices = as.integer(names(map(dataModel, 1))),
    Model = map(dataModel, 1) %>% unlist(use.names = F),
    correlation = map(dataModel, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > R2SOD,
    fittingModel = map(dataModel, 2)
  )
  #  nest(ChooseModel = -groupIndices)

  setDT(dataModel)

  processing[
    groupIndices %in% dataModel$groupIndices,
    aboveCorFit := rep(dataModel$aboveMinCor, each = nDilutions)
  ]
  # processing <- left_join(processList$processing, dataModel, by = c("groupIndices"), suffix = c(".x", ".y"))# |>
  # unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)


  # assert_that(nrow(processList$BestModel) == remainingCompound)

  # assert_that(nrow(processList$dataFittingModel) == remainingCompound * n_dilution )

  discardCompoundFitting <- n_distinct(dataModel |> filter(aboveMinCor == FALSE) %>% dplyr::select(groupIndices))
  remainingCompoundFitting <- n_distinct(dataModel |> filter(aboveMinCor == TRUE) %>% dplyr::select(groupIndices))

  message(discardCompoundFitting, " Compounds have a low R^2 (less than ", R2SOD, ") and undergo a second outlier detection.\n")

  # 2.) Second Outlier Detection

  if (discardCompoundFitting > 0) {
    if (nCORE > 1) plan(multisession, workers = nCORE)
    dataSOD <- my_fcn(
      xs = 1:n_distinct(processing |> filter(aboveCorFit == FALSE) |> dplyr::select(groupIndices)),
      inputData = processing |> filter(aboveCorFit == FALSE),
      func = outlierDetection,
      x = X,
      y = Intensity,
      numboutlier = nDilutions
    ) |>
      ldply(.id = NULL)

    plan(sequential)
    dataSOD <- setDT(dataSOD)[outlier %in% TRUE, Comment := str_replace(Comment, "outlier$", "outlierSOD")]
    # assert_that(n_distinct(dataSOD$groupIndices) == n_distinct(processing$groupIndices))
    processing[, Comment := as.character(Comment)][IDintern %in% dataSOD[dataSOD$outlier %in% TRUE, IDintern],
      ":="(Comment = dataSOD[dataSOD$outlier %in% TRUE, Comment],
        color = dataSOD[dataSOD$outlier %in% TRUE, color],
        pch = dataSOD[dataSOD$outlier %in% TRUE, pch],
        outlierSOD = dataSOD[dataSOD$outlier %in% TRUE, outlier]),
      on = "IDintern"
    ]


    message("For ", n_distinct(dataSOD |> dplyr::filter(outlier %in% TRUE) %>% dplyr::select(groupIndices)), " Compounds an Outlier were found\n")


    # assert_that(n_distinct(dataOutSec$groupIndices) == discardCompoundFitting)


    # 3. check length of consecutive points
    message("check length of consecutive points after second outlier detection")

    dataConsSOD <- processing[groupIndices %in% processing[outlierSOD %in% TRUE, groupIndices] & color %in% "black", .N, by = .(groupIndices)]
    dataConsSOD[, enoughPeaks := N >= MIN_FEATURE][enoughPeaks %in% F, Comment := "notEnoughPeaksSOD"]

    # assert_that(n_distinct(dataCons$groupIndices) == n_distinct(processing$groupIndices))

    processing[
      groupIndices %in% dataConsSOD[enoughPeaks %in% FALSE, groupIndices],
      enoughPeaks := rep(FALSE, each = nDilutions)
    ]

    dataConsSOD <- setDT(dataConsSOD)[, Comment := str_replace(Comment, "notEnoughPeaks", "notEnoughPeaksSOD")]

    processing[, Comment := as.character(Comment)][groupIndices %in% dataConsSOD[enoughPeaks %in% FALSE, groupIndices] & (!Comment %in% c(NA, NULL, "", " ")), Comment := paste(Comment, dataConsSOD[enoughPeaks %in% FALSE, Comment], sep = "_")]
    processing[, Comment := as.character(Comment)][groupIndices %in% dataConsSOD[enoughPeaks %in% FALSE, groupIndices] & (Comment %in% c(NA, NULL, "", " ")), Comment := "notEnoughPeaksSOD"]

    discardCompoundSOD <- n_distinct(dataConsSOD %>% filter(enoughPeaks == FALSE) %>% dplyr::select(groupIndices))
    remainingCompoundSOD <- n_distinct(dataConsSOD %>% filter(enoughPeaks != FALSE) %>% dplyr::select(groupIndices))
    message(discardCompound, " features had less than ", MIN_FEATURE, " Peaks and were discarded.\n--------------------------------------------------------\n")
    message(paste0("Data set with ", remainingCompound + remainingCompoundSOD, " Features.\n--------------------------------------------------------\n"))

    # 4.) choose model second time


    message("Fitting linear, logistic and quadratic regression after second outlierDetection\n")
    if (n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices]) > 0) {
      if (nCORE > 1) plan(multisession, workers = nCORE)
      dataModelSOD <- my_fcn2(
        xs = 1:n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices]),
        inputData = processing[groupIndices %in% processing[outlierSOD %in% TRUE & enoughPeaks %in% TRUE, groupIndices] & color %in% "black"],
        func = chooseModel,
        x = X,
        y = Intensity
      ) %>% unlist(recursive = F)


      dataModelSOD <- tibble(
        groupIndices = as.integer(names(map(dataModelSOD, 1))),
        Model = map(dataModelSOD, 1) %>% unlist(use.names = F),
        correlation = map(dataModelSOD, 3) %>% unlist(use.names = F),
        aboveMinCor = correlation > R2SOD,
        fittingModel = map(dataModelSOD, 2)
      )

      dataModelSOD <- setDT(dataModelSOD)

      processing[
        groupIndices %in% dataModelSOD[aboveMinCor %in% TRUE, groupIndices],
        aboveCorFit := rep(dataModelSOD[aboveMinCor %in% TRUE, aboveMinCor], each = nDilutions)
      ]
    } else {
      dataModelSOD <- NULL
    }

    discardCompoundFitting <- n_distinct(processing |> filter(aboveCorFit == FALSE) %>% dplyr::select(groupIndices))
    savedCompounds <- n_distinct(processing |> filter(aboveCorFit == TRUE, enoughPeaks %in% TRUE) %>% dplyr::select(groupIndices))
    message(savedCompounds, " of ", discardCompoundFitting, " Features have a R^2 above ", R2SOD, " after second outlier detection\n--------------------------------------------------------\n")

    # combine data for both outlier detections
    dataModelCombined <- copy(dataModel)
    if (n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices]) > 0) {
      dataModelCombined <- dataModelCombined[
        groupIndices %in% dataModelSOD[aboveMinCor %in% TRUE, groupIndices],
        ":="(Model = dataModelSOD[aboveMinCor %in% TRUE, Model],
          correlation = dataModelSOD[aboveMinCor %in% TRUE, correlation],
          aboveMinCor = dataModelSOD[aboveMinCor %in% TRUE, aboveMinCor],
          fittingModel = dataModelSOD[aboveMinCor %in% TRUE, fittingModel])
      ]
    }
  } else {
    processing$outlierSOD <- NA
    dataSOD <- NULL
    dataConsSOD <- NULL
    dataModelSOD <- NULL
    dataModelCombined <- copy(dataModel)
  }
  # assert_that(n_distinct(processList$processing |> filter(aboveMinCor %in% TRUE) |> dplyr::select(groupIndices)) == remainingCompoundFitting + savedCompounds )

  message("Data set with ", n_distinct(processing |> filter(aboveCorFit %in% TRUE) |> dplyr::select(groupIndices)), " Compounds with a R^2 above ", R2SOD, "\n--------------------------------------------------------\n")

  ## find linear Range
  message("Determining linear range\n--------------------------------------------------------\n")

  processing[groupIndices %in% dataModelCombined$groupIndices, fittingModel := rep(dataModelCombined$fittingModel, each = nDilutions)]

  if (nCORE > 1) plan(multisession, workers = nCORE)
  dataLinearRange <- my_fcn(
    xs = 1:n_distinct(processing |> filter(aboveCorFit %in% TRUE) |> dplyr::select(groupIndices)),
    inputData = processing |> filter(aboveCorFit %in% TRUE, color %in% "black"),
    func = findLinearRange,
    x = X,
    y = Intensity,
    modelObject = "fittingModel",
    MIN_FEATURE = MIN_FEATURE,
    res = res
  ) |>
    ldply(.id = NULL)

  plan(sequential)




  processinglinear <- left_join(processing[IDintern %in% dataLinearRange$IDintern, -c("pch", "color", "fittingModel")], dataLinearRange, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x, Comment.y), remove = TRUE, na.rm = TRUE)
  processing <- processing |> dplyr::select(-fittingModel)
  processing <- full_join(processing[!IDintern %in% processinglinear$IDintern], processinglinear, by = colnames(processing))


  # assert_that(n_distinct(processList$processing$groupIndices) == rawCompounds)

  message("For ", n_distinct(processing %>% filter(IslinearRange == TRUE & enoughPointsWithinLinearRange == TRUE) %>% dplyr::select(groupIndices)), " Features a linear Range with a minimum of ", MIN_FEATURE, " Points were found.\n--------------------------------------------------------\n")

  #
  # <!-- plotSignals(DAT = processList$Preprocessed$dataLinearRange %>% filter(ID %in% metabolitesRandomSample), x = "DilutionPoint", y = "IntensityNorm") -->
  #
  message("summarize all informations in one list\n--------------------------------------------------------\n")

  processing <- full_join(processing, dataPrep[!groupIndices %in% dropNA$groupIndices, -c("pch", "color", "N", "Comment", "enoughPeaks")],
    by = c("groupIndices", "IDintern", "Intensity", "X", "IntensityNorm", "DilutionPoint", "IntensityLog", "ConcentrationLog")
  )


  processing <- full_join(x = processing, y = dataPrep[groupIndices %in% dropNA$groupIndices, -c("N")], by = names("y"))


  setorder(processing, DilutionPoint)
  dataSummary <- getSummaryList(processing)

  processing <- dplyr::full_join(
    x = dataOrigin, y = processing,
    by = setNames(
      c("ID", "IDintern", "REPLICATE", "X", "Intensity"),
      c(COLNAMES[["ID"]], "IDintern", COLNAMES[["REPLICATE"]], COLNAMES[["X"]], COLNAMES[["Intensity"]])
    ),
    suffix = c(".raw", ".processed")
  )


  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  processList <- list(
    "dataOrigin_1" = dataOrigin,
    "dataPrep_2" = dataPrep,
    "dataTrim_4" = dataTrim,
    "dataCons_5" = dataCons,
    "dataModel_6" = dataModel,

    # "dataConsSOD_8" = dataConsSOD,
    # "dataModelSOD_9" = dataModelSOD,
    "dataLinearRange_10" = dataLinearRange,
    "summaryFFDS" = processing,
    "summaryFDS" = dataSummary,
    "Parameters" = list(
      data = DAT,
      COLNAMES = COLNAMES,
      nCORE = nCORE,
      calcFor = calcFor,
      SRES = SRES,
      R2FOD = R2FOD,
      R2SOD = R2SOD,
      MIN_FEATURE = MIN_FEATURE,
      X = X,
      Intensity = Intensity
    )
  )


  if (!empty(dataFOD)) processList$"dataFOD_3" <- dataFOD
  if (!is.null(dataSOD) | !empty(dataSOD)) processList$"dataSOD_7" <- dataSOD
  if (!is.null(dataConsSOD) | !empty(dataConsSOD)) processList$"dataConsSOD_8" <- dataConsSOD
  if (!is.null(dataModelSOD) | !empty(dataModelSOD)) processList$"dataModelSOD_9" <- dataModelSOD


  return(processList)
}

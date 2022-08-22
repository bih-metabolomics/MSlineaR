#' Title
#'
#' @param columnNames
#' @param dat
#' @param OutlierResPre
#' @param minCorFittingModelFOD
#' @param minCorFittingModel
#' @param minConsecutives
#'
#' @return
#' @export
#'
#' @examples
AssessLinearity <- function(columnNames,
                            dat,
                            nCore = 4,
                            OutlierResPre = 2,
                            minCorFittingModelFOD = 0.99,
                            minCorFittingModel = 0.99,
                            minConsecutives=5,
                            calcFor = c("Replicates", "Median", "both"),
                            Concentration = "Concentration",
                            Intensity = "Intensity"

){
  #progressbar
  pboptions(type = "timer", char = "[=-]", style = 5)
  #pb <- timerProgressBar(min = 0, max = grpn, width = 50, char = "[=-]", style = 5)

  progressr::handlers(global = TRUE)
  progressr::handlers(
    progressr::handler_progress(
      format   = ":current/:total [:bar] :percent in :elapsed ETA: :eta",
      width    = 60,
      complete = "+"
    ))

  #if(nCore > 1){


    #handlers(global = TRUE)
    my_fcn <- function(xs, func, inputData, ...) {

      p <- progressr::progressor(along = xs)
      y <- furrr::future_map(xs, function(x) {
        p(sprintf("x=%g", x))
        func(setDT(inputData)[groupIndices %in% unique(groupIndices)[x]], ...)
        #groupInd <- unique(inputData$groupIndices)[x]
        #inputData = filter(inputData, groupIndices %in% groupInd)
        #inputData[groupIndices %in% x, func(.SD, ...)]
        #func(inputData, ...)
      }, .progress = F)

    }

  # } else{
  my_fcn2 <- function(xs, func, inputData, ...) {
    p <- progressr::progressor(along = xs)
    y <- map(xs, function(x) {
      p(sprintf("x=%g", x))
      #inputData[groupIndices %in% x, func(.SD, ...)]
      func(setDT(inputData)[groupIndices %in% unique(groupIndices)[x]], ...)
      #groupInd <- unique(inputData$groupIndices)[x]
      #inputData = filter(inputData, groupIndices %in% groupInd)
      #func(inputData, ...)
    })
  }
  # }


  # Args: column names,file, residuals for outlier detection, minimal consecutive points in linear range, minimal fitting correlation
  dataOrigin <- data.table(dat)

  if(!"Replicate" %in% names(columnNames)){
    columnNames["Replicate"] = "Replicate"
    dataOrigin[, Replicate:= "1"]
    }

  setorderv(data.table(dataOrigin), c(columnNames[["ID"]], columnNames[["Replicate"]], columnNames[["Concentration"]]))
  dataOrigin[,IDintern:= paste0("s",1:.N)]

  #rm(dat)

  assert_that(all(columnNames %in% colnames(dataOrigin)), msg = "missing columns")

  #rename
  processing <- copy(dataOrigin)
  setnames(x = processing, old = columnNames[sort(names(columnNames))], new = c("Concentration", "ID", "Intensity", "mz", "Replicate", "rt"))
  processing <- processing[ , `:=`(Replicate=as.character(Replicate))][ , .(IDintern, ID, Replicate, mz, rt, Concentration, Intensity)]

  nCompounds = uniqueN(processing, by = c("ID"))
  nReplicates = uniqueN(processing, by = c("Replicate"))
  nDilutions = uniqueN(processing, by = c("Concentration"))
  nFeatures = uniqueN(processing, by = c("ID", "Replicate"))
  nPeaks = uniqueN(processing, by = c("IDintern"))

  message("Data set with ", nCompounds  ," molecular Features, ",
          nReplicates," Replicate(s) and ",
          nDilutions, " Concentration / Dilutions -> ",
          format.default(nFeatures, big.mark = ",", scientific = F), " Feature dilution series.",
          "\n--------------------------------------------------------\n")


  if(calcFor %in% c("both", "Median") & nReplicates >=2){

    message("calculating median and RSD of  ",nReplicates ," replicates per Compound\n--------------------------------------------------------\n")
    grpn = uniqueN(processing, by = c("ID", "Concentration"))
    combinedBatches <- copy(processing)
    combinedBatches[, IDintern := NULL][, `:=`(RSD =  sd(Intensity, na.rm = T)/mean(Intensity, na.rm = T)*100,
                                               Intensity = median(Intensity, na.rm = T),
                                               Replicate = "all"), by = .(ID, Concentration)]
    combinedBatches <-unique(combinedBatches)
    combinedBatches[,IDintern:= paste0("c",1:.N)]
    if(calcFor == "both"){
      processing <- dplyr::full_join(x = processing, y =  combinedBatches, by = intersect(names(combinedBatches), names(processing)))
    } else{
      processing <- combinedBatches
    }
rm("combinedBatches")

  }
# drop Compounds with values less than minConsecutives
  processing <- dplyr::left_join(processing[!is.na(Intensity), .N, by=.(ID, Replicate)][ N >= minConsecutives][ , N:= NULL ], processing , by = c("ID", "Replicate"))

  nCompoundsNew <- uniqueN(processing, by = c("ID"))
  nReplicatesNew <- uniqueN(processing, by = c("Replicate"))
  nFeaturesNew <- uniqueN(processing, by = c("ID", "Replicate"))

  message("Removed ",nCompounds - nCompoundsNew," molecular Features with less than ", minConsecutives, " points.
          New Data set with ", nCompoundsNew ," molecular Features and ",nReplicatesNew , " Replicates -> ",nFeaturesNew," Feature dilution series \n--------------------------------------------------------\n")

  ## normalizing, centralizing
  #message("normalising data\n--------------------------------------------------------\n")

  setorderv(processing, c("ID", "Replicate", "Concentration" ))

  processing[ , ':=' (Comment = NA, pch = fcase(!is.na(all_of(Intensity)), 19), color = fcase(!is.na(Intensity), "black"))]
  processing[ , ':=' (IntensityNorm = Intensity/max(Intensity, na.rm = T)*100,
                      DilutionPoint = 1:.N,
                      groupIndices = .GRP) ,by = c("ID", "Replicate")]


  dataPrep <- processing

  processing <- processing[, .(groupIndices, IDintern, Intensity, Concentration, IntensityNorm, DilutionPoint, Comment, pch, color)]

  #assert_that(nrow(processList$processing) == rawCompounds*n_dilution*nReplications)

  ##for testing
  #processList$processing <- processList$processing |> filter(groupIndices %in% 1:50)


  ## outlier

  message("First Outlier Detection\n")

  #pblapply(1:uniqueN(processing, by = "groupIndices"), function(i) processing[groupIndices %in% i, outlierDetection(.SD)]) |> ldply()

  if(nCore > 1) plan(multisession, workers = nCore)
   dataFOD <- my_fcn(xs = 1 : n_distinct(processing$groupIndices),
                     inputData = processing,
                     x = Concentration,
                     y = Intensity,
                     func = outlierDetection,
                     numboutlier = 1,
                     res = OutlierResPre,
                     threshCor = minCorFittingModelFOD) |>
     ldply(.id = NULL)

   plan(sequential)
   setDT(dataFOD)[outlier %in% TRUE, Comment := str_replace(Comment, "outlier$", "outlierFOD")]
   #assert_that(n_distinct(dataFOD$groupIndices) == n_distinct(processing$groupIndices))
   processing[, Comment := as.character(Comment)][IDintern %in% dataFOD$IDintern, ':=' (Comment = dataFOD$Comment,
                                                    color = dataFOD$color,
                                                    pch = dataFOD$pch,
                                                    outlierFOD = dataFOD$outlier), on = "IDintern"]


  message("For ",n_distinct(processing |> dplyr::filter(outlierFOD %in% TRUE) %>% dplyr::select(groupIndices))," Compounds an Outlier were found\n")

  ## trim
  message("Trim data: first Dilution should have the smallest Intensity and last point should have the biggest.")

  if(nCore > 1) plan(multisession, workers = nCore)
  dataTrim <- my_fcn(xs = 1 : n_distinct(processing$groupIndices),
                    inputData = processing,
                    func = trimEnds,
                    x = Concentration,
                    y = Intensity) |>
    ldply(.id = NULL)

  plan(sequential)

  #assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processing$groupIndices))

  processing[, Comment := as.character(Comment)][IDintern %in% dataTrim$IDintern, ':=' (Comment = dataTrim$Comment,
                                                                                       color = dataTrim$color,
                                                                                       pch = dataTrim$pch,
                                                                                       outlierFOD = dataTrim$outlierFOD,
                                                                                       trim = dataTrim$trim), on = "IDintern"]

  # check length of consecutive points
  message("check length of consecutive points")
TrimFeatures = n_distinct(dataTrim |> dplyr::filter(!color %in% "black") %>% dplyr::select(groupIndices))
message(TrimFeatures, "Features were removed after trimming the data")

  dataCons <- processing[color == "black", .N, by=.(groupIndices)]
  dataCons[, enoughPeaks:= N >= minConsecutives][enoughPeaks %in% F, Comment := "notEnoughPeaks"]

   #assert_that(n_distinct(dataCons$groupIndices) == n_distinct(processing$groupIndices))

  processing[groupIndices %in% dataCons$groupIndices,
               enoughPeaks := rep(dataCons$enoughPeaks, each = 9)]

  processing[ , Comment := as.character(Comment)][groupIndices %in% dataCons$groupIndices & (!Comment %in% c(NA, NULL, "", " ")), Comment := paste(Comment,dataCons$Comment, sep = "_")]
  processing[ , Comment := as.character(Comment)][groupIndices %in% dataCons$groupIndices & enoughPeaks %in% FALSE & (Comment %in% c(NA, NULL, "", " ")), Comment := "notEnoughPeaks"]




  discardCompound <- n_distinct(processing %>%  filter(enoughPeaks == FALSE, !color %in% "black") %>% dplyr::select(groupIndices))
  remainingCompound <- n_distinct(processing %>%  filter(enoughPeaks == TRUE, color %in% "black") %>% dplyr::select(groupIndices))
  message(discardCompound, " Compounds had less than ", minConsecutives, " Peaks and were discarded.\n--------------------------------------------------------\n")
  message("Data set with ", remainingCompound, " Compounds.\n--------------------------------------------------------\n")

  ##Fitting

  #1.) choose Model according to correlation threshold

  message("Fitting linear, logistic and quadratic regression\n")

  if(nCore > 1) plan(multisession, workers = nCore)
  dataModel <- my_fcn(xs = 1 : n_distinct(processing |> filter(enoughPeaks %in% TRUE, color %in% "black") |> dplyr::select(groupIndices)),
                     inputData = processing[enoughPeaks %in% TRUE &  color %in% "black"],
                     func = chooseModel,
                     x = Concentration,
                     y = Intensity

  ) %>% unlist(recursive = F)

  plan(sequential)

#gc()

  dataModel <- tibble(
    groupIndices = as.integer(names(map(dataModel, 1))),
    Model = map(dataModel, 1) %>% unlist(use.names = F),
    correlation = map(dataModel, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > minCorFittingModel,
    fittingModel = map(dataModel, 2))
  #  nest(ChooseModel = -groupIndices)

  setDT(dataModel)

  processing[groupIndices %in% dataModel$groupIndices,
             aboveCorFit := rep(dataModel$aboveMinCor, each = 9)]
  #processing <- left_join(processList$processing, dataModel, by = c("groupIndices"), suffix = c(".x", ".y"))# |>
    #unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)


  #assert_that(nrow(processList$BestModel) == remainingCompound)

  #assert_that(nrow(processList$dataFittingModel) == remainingCompound * n_dilution )

  discardCompoundFitting <- n_distinct(dataModel |>  filter(aboveMinCor == FALSE) %>% dplyr::select(groupIndices))
  remainingCompoundFitting <- n_distinct(dataModel |>  filter(aboveMinCor == TRUE) %>% dplyr::select(groupIndices))

  message(discardCompoundFitting, " Compounds have a low R^2 (less than " ,minCorFittingModel, ") and undergo a second outlier detection.\n" )

  #2.) Second Outlier Detection


  if(nCore > 1) plan(multisession, workers = nCore)
  dataSOD <- my_fcn(xs = 1 : n_distinct(processing|> filter(aboveCorFit == FALSE) |> dplyr::select(groupIndices)),
                    inputData = processing |> filter(aboveCorFit == FALSE),
                    func = outlierDetection,
                    x = Concentration,
                    y = Intensity,
                    numboutlier = nDilutions
                    ) |>
    ldply(.id = NULL)

  plan(sequential)
  dataSOD <- setDT(dataSOD)[outlier %in% TRUE, Comment := str_replace(Comment, "outlier$", "outlierSOD")]
  #assert_that(n_distinct(dataSOD$groupIndices) == n_distinct(processing$groupIndices))
  processing[, Comment := as.character(Comment)][IDintern %in% dataSOD[dataSOD$outlier %in% TRUE, IDintern],
                                                 ':=' (Comment = dataSOD[dataSOD$outlier %in% TRUE,Comment],
                                                       color = dataSOD[dataSOD$outlier %in% TRUE,color],
                                                       pch = dataSOD[dataSOD$outlier %in% TRUE,pch],
                                                       outlierSOD = dataSOD[dataSOD$outlier %in% TRUE,outlier]), on = "IDintern"]


  message("For ",n_distinct(dataSOD |> dplyr::filter(outlier %in% TRUE) %>% dplyr::select(groupIndices))," Compounds an Outlier were found\n")


  # assert_that(n_distinct(dataOutSec$groupIndices) == discardCompoundFitting)


  # 3. check length of consecutive points
  message("check length of consecutive points after second outlier detection")

  dataConsSOD <- processing[groupIndices %in% processing[outlierSOD %in% TRUE, groupIndices] & color %in% "black", .N, by=.(groupIndices)]
  dataConsSOD[, enoughPeaks:= N >= minConsecutives][enoughPeaks %in% F, Comment := "notEnoughPeaksSOD"]

  #assert_that(n_distinct(dataCons$groupIndices) == n_distinct(processing$groupIndices))

  processing[groupIndices %in% dataConsSOD[enoughPeaks %in% FALSE, groupIndices],
             enoughPeaks := rep(FALSE, each = 9)]

  dataConsSOD <- setDT(dataConsSOD)[, Comment := str_replace(Comment, "notEnoughPeaks", "notEnoughPeaksSOD")]

  processing[ , Comment := as.character(Comment)][groupIndices %in% dataConsSOD[enoughPeaks %in% FALSE, groupIndices] & (!Comment %in% c(NA, NULL, "", " ")), Comment := paste(Comment,dataConsSOD[enoughPeaks %in% FALSE, Comment], sep = "_")]
  processing[ , Comment := as.character(Comment)][groupIndices %in% dataConsSOD[enoughPeaks %in% FALSE, groupIndices] & (Comment %in% c(NA, NULL, "", " ")), Comment := "notEnoughPeaksSOD"]

  discardCompoundSOD <- n_distinct(dataConsSOD %>%  filter(enoughPeaks == FALSE) %>% dplyr::select(groupIndices))
  remainingCompoundSOD <- n_distinct(dataConsSOD %>%  filter(enoughPeaks != FALSE) %>% dplyr::select(groupIndices))
  message(discardCompound, " features had less than ", minConsecutives, " Peaks and were discarded.\n--------------------------------------------------------\n")
  message(paste0("Data set with ", remainingCompound + remainingCompoundSOD, " Features.\n--------------------------------------------------------\n"))

  # 4.) choose model second time


  message("Fitting linear, logistic and quadratic regression after second outlierDetection\n")
if(n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices])> 0){

  if(nCore > 1) plan(multisession, workers = nCore)
  dataModelSOD<- my_fcn2(xs = 1 : n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices]),
                       inputData = processing[groupIndices %in% processing[outlierSOD %in% TRUE & enoughPeaks %in% TRUE, groupIndices] & color %in% "black"],
                       func = chooseModel,
                       x = Concentration,
                       y = Intensity

  ) %>% unlist(recursive = F)


  dataModelSOD<- tibble(
    groupIndices = as.integer(names(map(dataModelSOD, 1))),
    Model = map(dataModelSOD, 1) %>% unlist(use.names = F),
    correlation = map(dataModelSOD, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > minCorFittingModel,
    fittingModel = map(dataModelSOD, 2))

  dataModelSOD <- setDT(dataModelSOD)

  processing[groupIndices %in% dataModelSOD[aboveMinCor %in% TRUE, groupIndices],
             aboveCorFit := rep(dataModelSOD[aboveMinCor %in% TRUE, aboveMinCor], each = 9)]

}

  discardCompoundFitting <- n_distinct(processing |>  filter(aboveCorFit == FALSE) %>% dplyr::select(groupIndices))
  savedCompounds <- n_distinct(processing |> filter(aboveCorFit == TRUE, enoughPeaks %in% TRUE) %>% dplyr::select(groupIndices))
  message(savedCompounds," of ", discardCompoundFitting," Features have a R^2 above ",minCorFittingModel ," after second outlier detection\n--------------------------------------------------------\n")

  # combine data for both outlier detections
  dataModelCombined <- copy(dataModel)
  if(n_distinct(dataConsSOD[enoughPeaks %in% TRUE, groupIndices])> 0){

  dataModelCombined <- dataModelCombined[groupIndices %in% dataModelSOD[aboveMinCor %in% TRUE, groupIndices],
            ':=' (Model = dataModelSOD[aboveMinCor %in% TRUE, Model],
                  correlation = dataModelSOD[aboveMinCor %in% TRUE, correlation],
                  aboveMinCor = dataModelSOD[aboveMinCor %in% TRUE, aboveMinCor],
                  fittingModel = dataModelSOD[aboveMinCor %in% TRUE, fittingModel])]
  }

  #assert_that(n_distinct(processList$processing |> filter(aboveMinCor %in% TRUE) |> dplyr::select(groupIndices)) == remainingCompoundFitting + savedCompounds )

  message("Data set with ", n_distinct(processing |> filter(aboveCorFit %in% TRUE) |> dplyr::select(groupIndices)), " Compounds with a R^2 above ",minCorFittingModel,"\n--------------------------------------------------------\n")

  ## find linear Range
  message("Determining linear range\n--------------------------------------------------------\n")

  processing[groupIndices %in% dataModelCombined$groupIndices, fittingModel := rep(dataModelCombined$fittingModel, each = 9)]

  if(nCore > 1) plan(multisession, workers = nCore)
  dataLinearRange <- my_fcn(xs = 1 : n_distinct(processing|> filter(aboveCorFit %in% TRUE) |> dplyr::select(groupIndices)),
                    inputData = processing |> filter(aboveCorFit %in% TRUE, color %in% "black"),
                    func = findLinearRange,
                    x = Concentration,
                    y = Intensity,
                    modelObject = "fittingModel",
                    minConsecutives = minConsecutives
  ) |>
    ldply(.id = NULL)

  plan(sequential)




 processinglinear <- left_join(processing[IDintern %in% dataLinearRange$IDintern, -c("pch", "color", "fittingModel")], dataLinearRange, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)
processing <- processing |> dplyr::select(-fittingModel)
processing <- full_join(processing[!IDintern %in% processinglinear$IDintern], processinglinear, by = colnames(processing))


  #assert_that(n_distinct(processList$processing$groupIndices) == rawCompounds)

  message("For ", n_distinct(processing %>% filter(IslinearRange == TRUE & enoughPointsWithinLinearRange == TRUE) %>% dplyr::select(groupIndices)), " Features a linear Range with a minimum of " ,minConsecutives ," Points were found.\n--------------------------------------------------------\n")

  #
  # <!-- plotSignals(dat = processList$Preprocessed$dataLinearRange %>% filter(ID %in% metabolitesRandomSample), x = "DilutionPoint", y = "IntensityNorm") -->
  #
  message("summarize all informations in one list\n--------------------------------------------------------\n")

  processing <- full_join(processing, dataPrep |> dplyr::select(ID, Replicate, IDintern, mz, rt), by = "IDintern")
  setorder(processing, DilutionPoint)
  dataSummary <- getSummaryList(processing)
  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  processList <- list("dataOrigin_1" = dataOrigin,
                      "dataPrep_2" = dataPrep,
                      "dataFOD_3" = dataFOD,
                      "dataTrim_4" = dataTrim,
                      "dataCons_5" = dataCons,
                      "dataModel_6" = dataModel,
                      "dataSOD_7" = dataSOD,
                      "dataConsSOD_8" = dataConsSOD,
                      "dataModelSOD_9" = dataModelSOD,
                      "dataLinearRange_10" = dataLinearRange,
                      "summaryFFDS" = processing,
                      "summaryFDS" = dataSummary)
  return(processList)

}


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
                            minConsecutives=5
){
  #progressbar
  pboptions(type = "timer", char = "-")

  progressr::handlers(global = TRUE)
  progressr::handlers(
    progressr::handler_progress(
      format   = ":current/:total [:bar] :percent in :elapsed ETA: :eta",
      width    = 60,
      complete = "+"
    ))

  if(nCore > 1){


    #handlers(global = TRUE)
    my_fcn <- function(xs, func, inputData, ...) {

      p <- progressr::progressor(along = xs)
      y <- furrr::future_map(xs, function(x) {
        p(sprintf("x=%g", x))
        groupInd <- unique(inputData$groupIndices)[x]
        inputData = filter(inputData, groupIndices %in% groupInd)
        func(inputData, ...)
      }, .progress = F)

    }

  } else{
  my_fcn <- function(xs, func, inputData, ...) {
    p <- progressr::progressor(along = xs)
    y <- map(xs, function(x) {
      p(sprintf("x=%g", x))
      groupInd <- unique(inputData$groupIndices)[x]
      inputData = filter(inputData, groupIndices %in% groupInd)
      func(inputData, ...)
    })
  }
  }


  # Args: column names,file, residuals for outlier detection, minimal consecutive points in linear range, minimal fitting correlation
processList <- list( dataOrigin = tibble(dat) |>
                       dplyr::arrange(get(columnNames[["ID"]]), get(columnNames[["Replicate"]]), get(columnNames[["Concentration"]])  ) |>
                       dplyr::mutate(IDintern = paste0("o",1:nrow(dat))))


  assert_that(all(columnNames %in% colnames(dat)), msg = "missing columns")

  #rename

  processList$processing <- processList$dataOrigin %>%
    dplyr::rename(ID = columnNames[["ID"]],
                  Replicate = ifelse(test = !is.null(columnNames[["Replicate"]]),yes = columnNames[["Replicate"]], no = "Replicate"),
                  mz = columnNames[["MZ"]],
                  RT = columnNames[["RT"]],
                  Concentration = columnNames[["Concentration"]],
                  Intensity = columnNames[["Intensity"]]) %>%
    dplyr::mutate(Replicate = ifelse(is.na(Replicate), "o1", as.character(Replicate)),
                  RSD = NA) |>
    dplyr::select(IDintern, ID, Replicate, mz, RT, Concentration, Intensity, RSD)



  #filter(Replicate %in% 1) %>%
  #mutate(IDintern = 1:nrow(.))

  message("Data set with ", n_distinct(processList$processing$ID) ," Compounds, ",
             n_distinct(processList$processing$Replicate)," Replicate(s) and ",
             n_distinct(processList$processing$Concentration),
             " Concentration / Dilutions\n--------------------------------------------------------\n")


  if(!is.null(columnNames[["Replicate"]]) & n_distinct(processList$processing$Replicate) >=2){

    message("calculating median and RSD of  ",n_distinct(processList$processing$Replicate) ," replicates per Compound\n--------------------------------------------------------\n")


  combinedBatches <-  processList$processing %>%
    dplyr::select(-IDintern) |>
    dplyr::group_by(ID, Concentration) %>%
    dplyr::mutate(
      RSD = sd(Intensity, na.rm = T)/mean(Intensity, na.rm = T)*100,
      Intensity = median(Intensity, na.rm = T),
      Replicate = "all") %>%
    distinct() |>
    ungroup() |>
    dplyr::mutate(IDintern = paste0("c",1:(n_distinct(processList$processing$ID)* n_distinct(processList$processing$Concentration))))

  processList$processing <- dplyr::full_join(x = processList$processing, y =  combinedBatches, by = c("IDintern", "ID", "Replicate", "mz", "RT", "Concentration", "Intensity", "RSD"))


}
  processList$processing <- processList$processing |>
    dplyr:: group_by(ID, Replicate) |>
    dplyr:: mutate( check = (!is.na(Intensity))) |>
    dplyr:: filter(sum(check) >= minConsecutives) |>
    dplyr::select(-check)

  rawCompounds <- n_distinct(processList$processing |>
                               dplyr:: group_by(ID, Replicate) |>
                               n_groups())
  nReplications <- n_distinct(processList$processing$Replicate)
  message("Removed compounds with less than ", minConsecutives, " values.\n
           Data set with ", n_distinct(processList$processing$ID) ," Compounds\n--------------------------------------------------------\n")

  n_dilution = n_distinct(processList$processing$Concentration)
  ## normalizing, centralizing
  message("normalising data\n--------------------------------------------------------\n")
  dataPrep <-   processList$processing   |>  prepareData()
  processList$processing <- full_join(x = dataPrep, y = processList$processing, by = "IDintern") |>
    nest(plotRawData = c(pch, color)) |>
    group_by(ID, Replicate) |>
    dplyr::mutate(groupIndices = cur_group_id()) |> ungroup()
  #splitting per batch

  #processList$processing <- processList$processing |> group_by(groupIn) #|> group_split()

  #assert_that(nrow(processList$processing) == rawCompounds*n_dilution*nReplications)

  ##for testing
  processList$processing <- processList$processing |> filter(groupIndices %in% 1:50)


  ## outlier

  message("First Outlier Detection\n")

  if(nCore > 1) plan(multisession, workers = nCore)
   dataOut <- my_fcn(xs = 1 : n_distinct(processList$processing$groupIndices),
                     inputData = processList$processing |>  unnest(cols = c(plotRawData)) |>  ungroup() |> group_by(groupIndices),
                     func = outlierDetection,
                     numboutlier = 1,
                     res = OutlierResPre,
                     threshCor = minCorFittingModelFOD) |>
     ldply(.id = NULL) |>
     as_tibble() |>
     dplyr::rename(outlierFOD = outlier) |>
     nest(plotFirstOutlierDetection = c(pch, color, residuals, ModelFit, correlationModel))


  plan(sequential)

  assert_that(n_distinct(dataOut$groupIndices) == n_distinct(processList$processing$groupIndices))

  processList$processing <- left_join(processList$processing, dataOut, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)

  message("For ",n_distinct(processList$processing |> dplyr::filter(outlierFOD %in% TRUE) %>% dplyr::select(groupIndices))," Compounds an Outlier were found\n")

  ## trim
  message("Trim data: first Dilution should have the smallest Intensity and last point should have the biggest.")

  if(nCore > 1) plan(multisession, workers = nCore)
  dataTrim <- my_fcn(xs = 1 : n_distinct(processList$processing$groupIndices),
                    inputData = processList$processing |> unnest(cols = plotFirstOutlierDetection),
                    func = trimEnds) |>
    ldply(.id = NULL) |>
    as_tibble() |>
    nest(plotTrimData = c(pch, color))

  plan(sequential)

  assert_that(n_distinct(dataTrim$groupIndices) == n_distinct(processList$processing$groupIndices))

  processList$processing <- left_join(processList$processing, dataTrim, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)



  # check length of consecutive points
  message("check length of consecutive points")

  if(nCore > 1) plan(multisession, workers = nCore)
  dataCons <- my_fcn(xs = 1 : n_distinct(processList$processing$groupIndices),
                     inputData = processList$processing |> unnest(cols = plotTrimData),
                     func = consecutiveVali,
                     minConsecutives = minConsecutives
                     ) |>
    ldply(.id = NULL) |>
    as_tibble()

  plan(sequential)

  assert_that(n_distinct(dataCons$groupIndices) == n_distinct(processList$processing$groupIndices))

  processList$processing <- left_join(processList$processing, dataCons, by = c("groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)


  discardCompound <- n_distinct(processList$processing %>%  filter(enoughPoints == FALSE) %>% dplyr::select(groupIndices))
  remainingCompound <- n_distinct(processList$processing %>%  filter(enoughPoints != FALSE) %>% dplyr::select(groupIndices))
  message(discardCompound, " Compounds had less than ", minConsecutives, " consecutive points and were discarded.\n--------------------------------------------------------\n")
  message(paste0("Data set with ", remainingCompound, " Compounds.\n--------------------------------------------------------\n"))

  ##Fitting

  #1.) choose Model according to correlation threshold

  message("Fitting linear, logistic and quadratic regression\n")

  if(nCore > 1) plan(multisession, workers = nCore)
  dataModel <- my_fcn(xs = 1 : n_distinct(processList$processing |> unnest(cols = plotTrimData) |> filter(enoughPoints %in% TRUE, color %in% "black") |> dplyr::select(groupIndices)),
                     inputData = processList$processing |> unnest(cols = plotTrimData) |> filter(enoughPoints %in% TRUE, color %in% "black"),
                     func = chooseModel
  ) %>% unlist(recursive = F)

  plan(sequential)


  dataModel <- tibble(
    groupIndices = as.integer(names(map(dataModel, 1))),
    Model = map(dataModel, 1) %>% unlist(use.names = F),
    correlation = map(dataModel, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > minCorFittingModel,
    fittingModel = map(dataModel, 2)) |>
    nest(ChooseModel = -groupIndices)

  processList$processing <- left_join(processList$processing, dataModel, by = c("groupIndices"), suffix = c(".x", ".y"))# |>
    #unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)


  #assert_that(nrow(processList$BestModel) == remainingCompound)

  #assert_that(nrow(processList$dataFittingModel) == remainingCompound * n_dilution )

  discardCompoundFitting <- n_distinct(processList$processing %>% unnest(cols = ChooseModel) |>  filter(aboveMinCor == FALSE) %>% dplyr::select(groupIndices))
  remainingCompoundFitting <- n_distinct(processList$processing %>% unnest(cols = ChooseModel) |>  filter(aboveMinCor == TRUE) %>% dplyr::select(groupIndices))

  message(discardCompoundFitting, " Compounds have a low R^2 (less than " ,minCorFittingModel, ") and undergo a second outlier detection.\n" )

  #2.) Second Outlier Detection
  dataSOD <- processList$processing %>% unnest(cols = c(ChooseModel)) |>  filter(aboveMinCor == FALSE)
  assert_that(nrow(dataSOD) == discardCompoundFitting * n_dilution)

  if(nCore > 1) plan(multisession, workers = nCore)
  dataOutSec <- my_fcn(xs = 1 : n_distinct(dataSOD$groupIndices),
                    inputData = dataSOD |>  unnest(cols = c(plotTrimData)) |>  filter(color %in% "black") |> ungroup() |> group_by(groupIndices),
                    func = outlierDetection,
                    numboutlier = n_dilution) |>
    ldply(.id = NULL) |>
    as_tibble() |>
    dplyr::rename(outlierSOD = outlier) |>
    nest(plotSecondOutlierDetection = c(pch, color, residuals, ModelFit, correlationModel))


  plan(sequential)

  assert_that(n_distinct(dataOutSec$groupIndices) == discardCompoundFitting)

  processList$processing <- left_join(processList$processing, dataOutSec, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)


  # 3.) choose model second time

  dataFittingModelSOD <-  processList$processing %>%  unnest(cols = plotSecondOutlierDetection) |> filter(color %in% "black")

  #1.) choose Model

  if(nCore > 1) plan(multisession, workers = nCore)
  dataModel <- my_fcn(xs = 1 : n_distinct(dataFittingModelSOD$groupIndices),
                      inputData = dataFittingModelSOD,
                      func = chooseModel
  ) %>% unlist(recursive = F)

  plan(sequential)


  dataModel <- tibble(
    groupIndices = as.integer(names(map(dataModel, 1))),
    Model = map(dataModel, 1) %>% unlist(use.names = F),
    correlation = map(dataModel, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > minCorFittingModel,
    fittingModel = map(dataModel, 2)) |>
    nest(ChooseModelSOD = -groupIndices)

  processList$processing <- left_join(processList$processing, dataModel, by = c("groupIndices"), suffix = c(".x", ".y"))# |>
 #
  discardCompoundFitting <- n_distinct(processList$processing %>% unnest(cols = ChooseModelSOD) |>  filter(aboveMinCor == FALSE) %>% dplyr::select(groupIndices))

  savedCompounds <- n_distinct(processList$processing |> unnest(cols = ChooseModelSOD)  |>  filter(aboveMinCor == TRUE) %>% dplyr::select(groupIndices))
  message(savedCompounds," of ", discardCompoundFitting," Compounds have a R^2 above ",minCorFittingModel ," after second outlier detection\n--------------------------------------------------------\n")

  # combine data for both outlierdetections
  combinedOutlier <- processList$processing |>
    unnest(cols = c(ChooseModel, ChooseModelSOD), keep_empty = T, names_sep = "_") |>
    mutate(aboveMinCor = ChooseModel_aboveMinCor ==TRUE | ChooseModelSOD_aboveMinCor ==TRUE) |>
    dplyr::select(IDintern, groupIndices, aboveMinCor)

  processList$processing <- left_join(processList$processing, combinedOutlier, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y"))# |>

  combinedOutlier <- processList$processing |>
    unnest(cols = c(ChooseModel, ChooseModelSOD), keep_empty = T, names_sep = "_") |>
    dplyr::select(IDintern, groupIndices, ChooseModel_aboveMinCor, ChooseModelSOD_aboveMinCor,
                  aboveMinCor, ChooseModel_fittingModel, ChooseModelSOD_fittingModel)|>
    filter(aboveMinCor %in% TRUE) |>
    mutate(check = unlist(map(1:nrow(processList$processing |> filter(aboveMinCor %in% TRUE)), ~!is.null(unlist(ChooseModelSOD_fittingModel)[x])))) |>
    mutate(fittingModel = ifelse(test = check %in% TRUE,
                                  yes = ChooseModelSOD_fittingModel,
                                   no = ChooseModel_fittingModel)) |>
    dplyr::select(IDintern, groupIndices, fittingModel)

  processList$processing <- left_join(processList$processing, combinedOutlier, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y"))# |>


  combinedOutlierplotData <- processList$processing |>
    dplyr::select(IDintern, groupIndices, plotTrimData, plotSecondOutlierDetection, outlierSOD, outlierFOD) |>
    unnest(cols = c(plotTrimData, plotSecondOutlierDetection),names_sep = "_", keep_empty = T) |>
    mutate(color = ifelse(test = is.na(plotSecondOutlierDetection_color), yes = plotTrimData_color, no = plotSecondOutlierDetection_color),
           pch = ifelse(test = is.na(plotSecondOutlierDetection_pch), yes = plotTrimData_pch, no = plotSecondOutlierDetection_pch),
           outlier = ifelse(test = outlierSOD %in% TRUE, yes = outlierSOD, no = outlierFOD)) |>

    nest(plotDataOutlier = c(color, pch)) |>
    dplyr::select(IDintern, groupIndices, outlier, plotDataOutlier)

  processList$processing <- left_join(processList$processing, combinedOutlierplotData, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y"))# |>



  assert_that(n_distinct(processList$processing |> filter(aboveMinCor %in% TRUE) |> dplyr::select(groupIndices)) == remainingCompoundFitting + savedCompounds )

  message("Data set with ", n_distinct(processList$processing |> filter(aboveMinCor %in% TRUE) |> dplyr::select(groupIndices)), " Compounds with a R^2 above ",minCorFittingModel,"\n--------------------------------------------------------\n")

  ## find linear Range
  message("Determining linear range\n--------------------------------------------------------\n")

  if(nCore > 1) plan(multisession, workers = nCore)
  dataLinearRange <- my_fcn(xs = 1 : n_distinct(processList$processing |> filter(aboveMinCor %in% TRUE) |> dplyr::select(groupIndices)),
                       inputData = processList$processing |>
                         filter(aboveMinCor %in% TRUE)|>
                         dplyr::select(IDintern, groupIndices,IntensityNorm, DilutionPoint, fittingModel, plotDataOutlier) |>
                         unnest(cols = plotDataOutlier) |>
                         filter(color %in% "black"),
                       func = findLinearRange,
                       modelObject = "fittingModel",
                       minConsecutives = minConsecutives) |>
    ldply(.id = NULL) |>
    as_tibble()  |>
    nest(plotLinearRange = c(pch, color, modelFit, ablineFit, ablineRes))


  plan(sequential)

  processList$processing <- left_join(processList$processing, dataLinearRange, by = c("IDintern", "groupIndices"), suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)

  assert_that(n_distinct(processList$processing$groupIndices) == rawCompounds)

  message("For ", n_distinct(processList$processing %>% filter(IslinearRange == TRUE & enoughPointsWithinLinearRange == TRUE) %>% dplyr::select(groupIndices)), " Compounds were a linear Range with a minimum of " ,minConsecutives ," Points found.\n--------------------------------------------------------\n")

  #
  # <!-- plotSignals(dat = processList$Preprocessed$dataLinearRange %>% filter(ID %in% metabolitesRandomSample), x = "DilutionPoint", y = "IntensityNorm") -->
  #
  message("summarize all informations in one list\n--------------------------------------------------------\n")

  processList$Summary <- getSummaryList(processList$processing)
  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  return(processList)

}


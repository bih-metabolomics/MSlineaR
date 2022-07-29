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
                            OutlierResPre = 2,
                            minCorFittingModelFOD = 0.99,
                            minCorFittingModel = 0.99,
                            minConsecutives=5
){
  pboptions(type = "timer", char = "-")
  # Args: column names,file, residuals for outlier detection, minimal consecutive points in linear range, minimal fitting correlation

  processList <- list(dataOrigin = dat)

  assert_that(all(columnNames %in% colnames(dat)), msg = "missing columns")

  #rename

  processList$dataRaw <- processList$dataOrigin %>%
    dplyr::rename(ID = columnNames[["ID"]],
                  Replicate = ifelse(test = !is.null(columnNames[["Replicate"]]),yes = columnNames[["Replicate"]], no = "Replicate"),
                  mz = columnNames[["MZ"]],
                  RT = columnNames[["RT"]],
                  Concentration = columnNames[["Concentration"]],
                  Intensity = columnNames[["Intensity"]]) %>%
    dplyr::select(ID, Replicate, mz, RT, Concentration, Intensity)
  #filter(Replicate %in% 1) %>%
  #mutate(IDintern = 1:nrow(.))

  cat(paste0("Data set with ", n_distinct(processList$dataRaw$ID) ," Compounds, ",
             n_distinct(processList$dataRaw$Replicate)," Replicate(s) and ",
             n_distinct(processList$dataRaw$Concentration),
             " Concentration / Dilutions\n--------------------------------------------------------\n"))


  if(is.null(columnNames[["Replicate"]])){
    processList$dataRaw$Replicate = 1
  } else{
    cat(paste0("calculating median and RSD of  ",n_distinct(processList$dataRaw$Replicate) ," replicates per Compound\n--------------------------------------------------------\n"))
  }

  processList$dataRepMed <-  processList$dataRaw %>%
    dplyr::group_by(ID, Concentration) %>%
    dplyr::summarize(.groups = "keep", IntensityMedian = median(Intensity, na.rm = T),
                     RSD = sd(Intensity, na.rm = T)/mean(Intensity, na.rm = T)*100) %>%
    full_join(., processList$dataRaw %>% dplyr::select(-Replicate, -Intensity) %>% distinct(),
              by = c("ID", "Concentration")) %>%
    mutate(IDintern = 1:nrow(.))

  processList$dataRepMed <- processList$dataRepMed |>
    dplyr:: group_by(ID) |>
    dplyr:: mutate( check = (!is.na(IntensityMedian))) |>
    dplyr:: filter(sum(check) >= minConsecutives) |>
    dplyr::select(-check)

  cat(paste0("Removed compounds with less than ", minConsecutives, "values.\n
           Data set with ", n_distinct(processList$dataRepMed$ID) ," Compounds, ",
             n_distinct(processList$dataRepMed$Replicate)," Replicate(s) and ",
             n_distinct(processList$dataRepMed$Concentration),
             " Concentration / Dilutions\n--------------------------------------------------------\n"))


  n_dilution = n_distinct(processList$dataRepMed$Concentration)
  ## normalizing, centralizing
  cat(paste0("log transforming and normalising\n--------------------------------------------------------\n"))
  processList$dataPrep <- processList$dataRepMed  %>% prepareData(dat = .)

  assert_that(nrow(processList$dataPrep) == nrow(processList$dataRepMed))

  ## outlier

  cat(paste0("First Outlier Detection\n"))
  dataOut <- pblapply(unique(processList$dataPrep$ID),
                      function(x) outlierDetection(dat = processList$dataPrep %>% filter(ID %in% x),
                                                   numboutlier = 1, res = OutlierResPre,
                                                   threshCor = minCorFittingModelFOD)) %>%
    unlist(recursive = F)

  assert_that(length(dataOut) == length(unique(processList$dataPrep$ID)))
  # FOD = First Outlier detection
  processList$BestModelFOD <- map(dataOut, 1) %>% ldply(.id = NULL)
  assert_that(nrow(processList$BestModelFOD) == length(unique(processList$dataPrep$ID)))

  processList$OutlierFOD <- map(dataOut, 2) %>% ldply(.id = NULL)
  assert_that(nrow(processList$OutlierFOD) == length(unique(processList$dataPrep$ID[!is.na(processList$dataPrep$IntensityRaw)])))
  assert_that(all(duplicated(processList$OutlierFOD$ID)) == FALSE)

  processList$dataFOD <- map(dataOut, 3) %>% ldply(.id = NULL)
  assert_that(nrow(processList$dataFOD) == nrow(processList$dataPrep))

  cat(paste0("For ",n_distinct(processList$dataFOD %>%  filter(Outlier == TRUE) %>% dplyr::select(ID))," Compounds were an Outlier found\n"))

  ## trim
  processList$dataTrim <- pblapply(unique(processList$dataFOD$ID),
                                   function(x) trimEnds(dats = processList$dataFOD %>%
                                                          filter(ID %in% x))) %>%
    ldply(.id = NULL)
  assert_that(nrow(processList$dataTrim) == nrow(processList$dataPrep))

  # check length of consecutive points
  processList$dataTrim <- pblapply(unique(processList$dataTrim$ID), function(x) {consecutiveVali(dats = processList$dataTrim %>% filter(ID %in% x), minConsecutives = 5)}) %>% ldply
  assert_that(nrow(processList$dataTrim) == nrow(processList$dataPrep))

  discardCompound <- n_distinct(processList$dataTrim %>%  filter(enoughPoints == FALSE) %>% dplyr::select(ID))
  remainingCompound <- n_distinct(processList$dataTrim %>%  filter(enoughPoints != FALSE) %>% dplyr::select(ID))
  cat(discardCompound, "Compounds had less than", minConsecutives, "consecutive points and were discarded.\n--------------------------------------------------------\n")
  cat(paste0("Data set with ", remainingCompound, " Compounds.\n--------------------------------------------------------\n"))

  ##Fitting

  #1.) choose Model according to correlation threshold

  cat("Fitting linear, logistic and quadratic regression\n")

  processList$dataFittingModel <- processList$dataTrim %>% filter(enoughPoints %in% TRUE)

  dataModel <- pblapply(unique(processList$dataFittingModel$ID), function(x) chooseModel(dat = processList$dataFittingModel %>% filter(ID %in% x, color %in% "black" ))) %>% unlist(recursive = F)

  processList$BestModel <- tibble(
    ID = names(map(dataModel, 1)),
    Model = map(dataModel, 1) %>% unlist(use.names = F),
    correlation = map(dataModel, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > minCorFittingModel)
  assert_that(nrow(processList$BestModel) == remainingCompound)

  processList$fittingModel <- map(dataModel, 2)
  processList$dataFittingModel <- full_join(processList$dataFittingModel, processList$BestModel, by = "ID")
  assert_that(nrow(processList$dataFittingModel) == remainingCompound * n_dilution )

  discardCompoundFitting <- n_distinct(processList$dataFittingModel %>%  filter(aboveMinCor == FALSE) %>% dplyr::select(ID))
  remainingCompoundFitting <- n_distinct(processList$dataFittingModel %>%  filter(aboveMinCor == TRUE) %>% dplyr::select(ID))

  cat(discardCompoundFitting, "Compounds have a low R^2 (less than" ,minCorFittingModel, ") and undergo a second outlier detection.\n" )

  #2.) Second Outlier Detection
  dataSOD <- processList$dataFittingModel %>% filter(aboveMinCor == FALSE)
  assert_that(nrow(dataSOD) == discardCompoundFitting * n_dilution)

  dataOutSec <- pblapply(unique(dataSOD$ID), function(x) outlierDetection(dat = dataSOD %>% filter(ID %in% x), numboutlier = n_dilution)) %>% unlist(recursive = F)
  assert_that(length(dataOutSec) == discardCompoundFitting)


  processList$BestModelSOD <- map(dataOutSec, 1) %>% ldply(.id = NULL)
  processList$OutlierSOD <- map(dataOutSec, 2) %>% ldply(.id = NULL)
  processList$dataSOD <- map(dataOutSec, 3) %>% ldply(.id = NULL)



  # 3.) choose model second time

  processList$dataFittingModelSOD <- processList$dataSOD %>%  filter(color %in% "black", enoughPoints %in% TRUE)

  #1.) choose Model

  dataModel <- pblapply(unique(processList$dataFittingModelSOD$ID), function(x) chooseModel(dat = processList$dataFittingModelSOD %>% filter(ID %in% x))) %>% unlist(recursive = F)

  processList$BestModelFittingSOD <- tibble(
    ID = names(map(dataModel, 1)),
    Model = map(dataModel, 1) %>% unlist(use.names = F),
    correlation = map(dataModel, 3) %>% unlist(use.names = F),
    aboveMinCor = correlation > minCorFittingModel)
  processList$fittingModelSOD <- map(dataModel, 2)

  processList$dataFittingModelSOD <- full_join(processList$dataFittingModelSOD , processList$BestModelFittingSOD, by = "ID", suffix = c(".first", ".second") )

  savedCompounds <- n_distinct(processList$dataFittingModelSOD %>% filter(aboveMinCor.second == TRUE) %>% dplyr::select(ID))
  cat(savedCompounds," of ", discardCompoundFitting,"Compounds have a R^2 above",minCorFittingModel ,"after second outlier detection\n--------------------------------------------------------\n")

  processList$dataFittingModelAll <- bind_rows(processList$dataFittingModel %>%  filter(aboveMinCor == TRUE),
                                               processList$dataFittingModelSOD %>%  filter(aboveMinCor.second == TRUE) %>%
                                                 mutate(Model = Model.second,
                                                        correlation = correlation.second,
                                                        aboveMinCor = aboveMinCor.second) %>%
                                                 dplyr::select(-contains(c("first", "second"))), .id = NULL
  )
  processList$dataFittingModelAll <- bind_rows(processList$dataFittingModelAll, processList$dataSOD %>% filter(ID %in% processList$dataFittingModelAll$ID, !IDintern %in% processList$dataFittingModelAll$IDintern))

  assert_that(n_distinct(processList$dataFittingModelAll$ID) == remainingCompoundFitting + savedCompounds )

  cat(paste0("Data set with ", n_distinct(processList$dataFittingModelAll$ID), " Compounds with a R^2 above ",minCorFittingModel,"\n--------------------------------------------------------\n"))

  processList$fittingModelAll <- c(processList$fittingModel[names(processList$fittingModel) %in% unlist(distinct(processList$dataFittingModel %>% filter(aboveMinCor == TRUE) %>% dplyr::select(ID)))],
                                   processList$fittingModelSOD[names(processList$fittingModelSOD) %in% unlist(distinct(processList$dataFittingModelSOD %>% filter(aboveMinCor.second == TRUE) %>% dplyr::select(ID)))])

  ## find linear Range
  cat("Determining linear range\n--------------------------------------------------------\n")
  processList$dataLinearRange <- pblapply(unique(processList$dataFittingModelAll$ID), function(x) findLinearRange(dat = processList$dataFittingModelAll %>% filter(ID %in% x, color %in% "black"), modelObject = processList$fittingModelAll[[x]])) %>% ldply

  processList$dataLinearRange <- bind_rows(processList$dataLinearRange, processList$dataFittingModelAll %>% filter(ID %in% processList$dataLinearRange$ID, !IDintern %in% processList$dataLinearRange$IDintern))

  processList$dataLinearRange <- bind_rows(processList$dataLinearRange, processList$dataTrim %>% filter(!IDintern %in% processList$dataLinearRange$IDintern))

  cat("For ", n_distinct(processList$dataLinearRange %>% filter(linearRange == TRUE) %>% dplyr::select(ID)), "Compounds were a linear Range found.\n--------------------------------------------------------\n")

  #
  # <!-- plotSignals(dat = processList$Preprocessed$dataLinearRange %>% filter(ID %in% metabolitesRandomSample), x = "DilutionPoint", y = "IntensityNorm") -->
  #
  processList$Summary <- getSummaryList(processList)
  # <!-- processList$SummaryAll <- getAllList(processList) -->
  #

  return(processList)

}


#' Title
#'
#' @param completeList
#'
#' @return
#' @export
#'
#' @examples
getSummaryList <- function(completeList){
  #browser()
  dat <- completeList[["dataRepMed"]] %>%
    dplyr::ungroup() %>%
    dplyr::select(ID, mz, RT) %>%
    distinct(.)

  LR <- completeList[["dataLinearRange"]] %>%
    dplyr::select(ID, linearRangeStart, linearRangeEnd) %>%
    distinct(.) %>%
    mutate(LengthLinearRange = linearRangeEnd - linearRangeStart + 1) %>%
    drop_na(linearRangeStart)

  LR_limits <-  lapply(LR$ID, function(x) {completeList[["dataPrep"]] %>%
      filter(ID %in% x) %>%
      mutate(LR_ConcentrationStart = ConcentrationRaw[LR$linearRangeStart[LR$ID %in% x]],
             LR_ConcentrationEnd = ConcentrationRaw[LR$linearRangeEnd[LR$ID %in% x]],
             LR_IntensityStart = IntensityRaw[LR$linearRangeStart[LR$ID %in% x]],
             LR_IntensityEnd = IntensityRaw[LR$linearRangeEnd[LR$ID %in% x]]) %>%
      dplyr::select(ID, LR_ConcentrationStart, LR_ConcentrationEnd, LR_IntensityStart, LR_IntensityEnd) %>%
      distinct()
  }) %>% ldply


  MissingNr <- completeList[["dataRepMed"]] %>%
    dplyr::ungroup() %>%
    dplyr::select(ID, IntensityMedian) %>%
    group_by(ID) %>%
    dplyr::count(is.na(IntensityMedian)) %>%
    dplyr::rename(Nr_MissingValue = n) %>%
    filter(`is.na(IntensityMedian)` == TRUE) %>%
    dplyr::select(ID, Nr_MissingValue)

  Missing <- lapply(LR$ID, function(x){
    completeList[["dataRaw"]] %>%
      filter(ID %in% x) %>%
      filter(row_number() %in% LR$linearRangeStart[LR$ID %in% x]: LR$linearRangeEnd[LR$ID %in% x]) %>%
      summarise(ID, "MissingInLinearRange" = is.na(Intensity)) %>%
      distinct()
  }) %>% ldply

  OutlierNr <- completeList[["dataLinearRange"]] %>%
    filter(Outlier %in% TRUE) %>%
    group_by(ID) %>%
    dplyr::count(Outlier) %>%
    dplyr::rename(Nr_Outlier = n) %>%
    dplyr::select(ID, Nr_Outlier)

  Outlier <- lapply(LR$ID, function(x){
    completeList[["dataLinearRange"]] %>%
      filter(ID %in% x) %>%
      filter(row_number() %in% LR$linearRangeStart[LR$ID %in% x]: LR$linearRangeEnd[LR$ID %in% x]) %>%
      dplyr::summarise(ID, "OutlierInLinearRange" = any(Outlier %in% TRUE, na.rm = T)) %>%
      distinct()
  }) %>% ldply


  all <- join_all(list(dat, LR, LR_limits, MissingNr, Missing, OutlierNr, Outlier), by = "ID")
  #all$OutlierInLinearRange[!is.na(all$Nr_Outlier) & all$Nr_Outlier>0 & is.na(all$OutlierInLinearRange)] <- FALSE
  all$OutlierInLinearRange[is.na(all$Nr_Outlier)] <- NA
  all$MissingInLinearRange[is.na(all$Nr_MissingValue)] <- NA

  all <- all %>% distinct()


  return(all)

}

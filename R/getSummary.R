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
  dat <- completeList %>%
    dplyr::ungroup() %>%
    dplyr::select(groupIndices,ID, mz, RT, Replicate) %>%
    distinct(.)

  LR <- completeList %>%
    dplyr::select(groupIndices, linearRangeStart, linearRangeEnd, linearRange, enoughPointsWithinLinearRange) %>%
    distinct(.) %>%
    drop_na(linearRangeStart)

  LR_limits <-  map(LR$groupIndices, .f = function(x) {completeList %>%
      filter(groupIndices %in% x) %>%
      mutate(LR_ConcentrationStart = Concentration[LR$linearRangeStart[LR$groupIndices %in% x]],
             LR_ConcentrationEnd = Concentration[LR$linearRangeEnd[LR$groupIndices %in% x]],
             LR_IntensityStart = Intensity[LR$linearRangeStart[LR$groupIndices %in% x]],
             LR_IntensityEnd = Intensity[LR$linearRangeEnd[LR$groupIndices %in% x]]) %>%
      dplyr::select(groupIndices, LR_ConcentrationStart, LR_ConcentrationEnd, LR_IntensityStart, LR_IntensityEnd) %>%
      distinct()
  }) %>% ldply


  MissingNr <- completeList %>%
    dplyr::ungroup() %>%
    dplyr::select(groupIndices, Intensity) %>%
    group_by(groupIndices) %>%
    dplyr::count(is.na(Intensity)) %>%
    dplyr::rename(Nr_MissingValue = n) %>%
    filter(`is.na(Intensity)` == TRUE) %>%
    dplyr::select(groupIndices, Nr_MissingValue)

  Missing <- map(LR$groupIndices, function(x){
    completeList %>%
      filter(groupIndices %in% x) %>%
      filter(row_number() %in% LR$linearRangeStart[LR$groupIndices %in% x]: LR$linearRangeEnd[LR$groupIndices %in% x]) %>%
      summarise(groupIndices, "MissingInLinearRange" = is.na(Intensity)) %>%
      distinct()
  }) %>% ldply

  OutlierNr <- completeList %>%
    filter(outlier %in% TRUE) %>%
    group_by(groupIndices) %>%
    dplyr::count(outlier) %>%
    dplyr::rename(Nr_Outlier = n) %>%
    dplyr::select(groupIndices, Nr_Outlier)

  Outlier <- map(LR$groupIndices, function(x){
    completeList %>%
      filter(groupIndices %in% x) %>%
      filter(row_number() %in% LR$linearRangeStart[LR$groupIndices %in% x]: LR$linearRangeEnd[LR$groupIndices %in% x]) %>%
      dplyr::summarise(groupIndices, "OutlierInLinearRange" = any(outlier %in% TRUE, na.rm = T)) %>%
      distinct()
  }) %>% ldply


  all <- join_all(list(dat, LR, LR_limits, MissingNr, Missing, OutlierNr, Outlier), by = "groupIndices")
  #all$OutlierInLinearRange[!is.na(all$Nr_Outlier) & all$Nr_Outlier>0 & is.na(all$OutlierInLinearRange)] <- FALSE
  all$OutlierInLinearRange[is.na(all$Nr_Outlier)] <- NA
  all$MissingInLinearRange[is.na(all$Nr_MissingValue)] <- NA

  all <- all %>% distinct()


  return(all)

}

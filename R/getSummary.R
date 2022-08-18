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

  dat <- setorder(na.last = T, completeList[ , .(groupIndices,mz, rt, Replicate, Comment)], Comment) |> distinct(.keep_all = T, groupIndices,mz, rt, Replicate)

  LR <- unique(completeList[!is.na(linearRangeStart) , .(groupIndices, linearRangeStart, linearRangeEnd, linearRange, enoughPointsWithinLinearRange)])

  LR_limits <-  map(LR$groupIndices, .f = function(x) {completeList %>%
      filter(groupIndices %in% x) %>%
      mutate(LR_ConcentrationStart = Concentration[LR$linearRangeStart[LR$groupIndices %in% x]],
             LR_ConcentrationEnd = Concentration[LR$linearRangeEnd[LR$groupIndices %in% x]],
             LR_IntensityStart = Intensity[LR$linearRangeStart[LR$groupIndices %in% x]],
             LR_IntensityEnd = Intensity[LR$linearRangeEnd[LR$groupIndices %in% x]]) %>%
      dplyr::select(groupIndices, LR_ConcentrationStart, LR_ConcentrationEnd, LR_IntensityStart, LR_IntensityEnd) %>%
      distinct()
  }) %>% ldply


  MissingNr <- unique(completeList[, .(groupIndices, Intensity)][,Nr_MissingValue := sum(is.na(Intensity)),"groupIndices"][,-("Intensity")])

  Missing <- map(LR$groupIndices, function(x){
    completeList %>%
      filter(groupIndices %in% x) %>%
      filter(row_number() %in% LR$linearRangeStart[LR$groupIndices %in% x]: LR$linearRangeEnd[LR$groupIndices %in% x]) %>%
      summarise("groupIndices" = unique(groupIndices), "MissingInLinearRange" = any(is.na(Intensity))) %>%
      distinct()
  }) %>% ldply

  OutlierNr <- completeList[outlierFOD %in% TRUE | outlierSOD %in% TRUE ,
                            .(IDintern,groupIndices, outlierFOD, outlierSOD)
                            ][, Nr_Outlier := sum(outlierFOD %in% TRUE |
                                                    outlierSOD %in% TRUE),
                              groupIndices][,.(groupIndices, Nr_Outlier)]


  Outlier <- map(LR$groupIndices, function(x){
    completeList %>%
      filter(groupIndices %in% x) %>%
      filter(row_number() %in% LR$linearRangeStart[LR$groupIndices %in% x]: LR$linearRangeEnd[LR$groupIndices %in% x]) %>%
      dplyr::summarise(groupIndices, "OutlierInLinearRange" = any(outlierFOD %in% TRUE | outlierSOD %in% TRUE, na.rm = T)) %>%
      distinct()
  }) %>% ldply


  all <- join_all(list(dat, LR, LR_limits, MissingNr, Missing, OutlierNr, Outlier), by = "groupIndices")
  #all$OutlierInLinearRange[!is.na(all$Nr_Outlier) & all$Nr_Outlier>0 & is.na(all$OutlierInLinearRange)] <- FALSE
  all$OutlierInLinearRange[is.na(all$Nr_Outlier)] <- NA
  all$MissingInLinearRange[is.na(all$Nr_MissingValue)] <- NA
  all$Comment[all$MissingInLinearRange %in% TRUE] <- paste0(all$Comment[all$MissingInLinearRange %in% TRUE], "_withMissing")
  all$Comment[all$OutlierInLinearRange %in% TRUE] <- paste0(all$Comment[all$MissingInLinearRange %in% TRUE], "_withOutlier")

  all$filter <- case_when(
    all$linearRange >= 1 & all$MissingInLinearRange %in% c(NA,FALSE) & all$OutlierInLinearRange %in% c(NA,FALSE) ~ "linear",
    all$linearRange >= 1 & all$MissingInLinearRange %in% TRUE & all$OutlierInLinearRange %in% c(NA,FALSE) ~ "linearWithMissing",
    all$linearRange >= 1 & all$MissingInLinearRange %in% c(NA,FALSE) & all$OutlierInLinearRange %in% TRUE ~ "linearWithOutlier",
    all$linearRange >= 1 & all$MissingInLinearRange %in% TRUE & all$OutlierInLinearRange %in% TRUE ~ "linearWithMissingAndOutlier",
    TRUE ~ "notLinear"
  )

  all <- all %>% distinct()


  return(all)

}

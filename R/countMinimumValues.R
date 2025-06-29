#' Title
#'
#' @param DAT
#' @param MIN_FEATURE minimum number of features within a dilution or concentration series to be considered as linear range
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
#' @import data.table
countMinimumValue <- function(DAT, MIN_FEATURE = parent.frame()$MIN_FEATURE, step, y){

  stopifnot(exprs ={"Parameter MIN_FEATURE needs to be greater or equal than 3" = MIN_FEATURE >= 3})

  dat <- data.table::data.table(DAT)

  #N <- paste0("N_",step)
  #enoughPeaks <- paste("enoughPeaks_",step)
  dat <- unique(dat[, N := sum(!is.na(get(y))),groupIndices][ , enoughPeaks := ifelse( N >= MIN_FEATURE, TRUE, FALSE)#,
                                                     #Comment = ifelse( N >= MIN_FEATURE, "EnoughPeaks", "notEnoughPeaks")
                                                     ][
                                                       ,list(groupIndices,Feature_ID, Batch, N, enoughPeaks)])

  DAT$color[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]] <- "grey"
  DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]] <- unlist(apply(cbind(DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]], "notEnoughPeaks"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))
  DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% TRUE]] <- unlist(apply(cbind(DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% TRUE]], "EnoughPeaks"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))
  DAT[[y]][DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]] <- NA

  colnames(dat)[which(colnames(dat) %in% "N")] <- paste0("N_",step)
  colnames(dat)[which(colnames(dat) %in% "enoughPeaks")] <- paste0("enoughPeaks_",step)

    return(list(DAT, dat))
}

#' Title
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
checkLength <- function(step, processingGroup, processingFeature, Compounds,
                        Dilutions, Series, Signals,
                        MIN_FEATURE){
  peaksOld <- ifelse(step > 1, paste0("enoughPeaks_", step-1), "enoughPeaks")

  peaksNew <- paste0("enoughPeaks_", step)

  nCompoundsOld <- data.table::uniqueN(processingGroup[get(peaksOld) %in% TRUE], by = c("Feature_ID"))
  nReplicatesOld <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("Batch"))
  nSeriesOld <- data.table::uniqueN(processingGroup[get(peaksOld) %in% TRUE], by = c("groupIndices"))

  nCompoundsNew <- data.table::uniqueN(processingGroup[get(peaksNew) %in% TRUE], by = c("Feature_ID"))
  nReplicatesNew <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("Batch"))
  nSeriesNew <- data.table::uniqueN(processingGroup[get(peaksNew) %in% TRUE], by = c("groupIndices"))

  logr::put(paste0("Removed ", nSeriesOld - nSeriesNew, " ", Series, " with less than ", MIN_FEATURE, " Signals."))
  logr::put(paste0("Remaining Data set with ", nCompoundsNew, " ", Compounds," and ", nReplicatesNew, " Batch(es) -> ", nSeriesNew, " ", Series))

  stopifnot(exprs = {
    "all Compounds were removed" = nCompoundsNew  > 0
    "all Dilution/Concentration-Series were removed" = nSeriesNew > 0
  })

  cutoff <- processingFeature[groupIndices %in% processingGroup[get(peaksNew) %in% FALSE, groupIndices]]
  processingFeature <- processingFeature[groupIndices %in% processingGroup[get(peaksNew) %in% TRUE, groupIndices]]

  return(list(processingFeature, cutoff))


}

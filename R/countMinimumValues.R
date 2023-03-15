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
countMinimumValue <- function(DAT, MIN_FEATURE = 3, step, y){

  stopifnot(exprs ={"Parameter MIN_FEATURE needs to be greater or equal than 3" = MIN_FEATURE >= 3})

  dat <- data.table::data.table(DAT)

  #N <- paste0("N_",step)
  #enoughPeaks <- paste("enoughPeaks_",step)
  dat <- unique(dat[, N := sum(color %in% "black"),groupIndices][ , enoughPeaks := ifelse( N >= MIN_FEATURE, TRUE, FALSE)#,
                                                     #Comment = ifelse( N >= MIN_FEATURE, "EnoughPeaks", "notEnoughPeaks")
                                                     ][
                                                       ,list(groupIndices,ID, Batch,  N, enoughPeaks)])

  DAT$color[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]] <- "grey"
  DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]] <- unlist(apply(cbind(DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]], "notEnoughPeaks"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))
  DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% TRUE]] <- unlist(apply(cbind(DAT$Comment[DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% TRUE]], "EnoughPeaks"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))
  DAT[[y]][DAT$groupIndices %in% dat$groupIndices[dat$enoughPeaks %in% FALSE]] <- NA

  colnames(dat)[which(colnames(dat) %in% "N")] <- paste0("N_",step)
  colnames(dat)[which(colnames(dat) %in% "enoughPeaks")] <- paste0("enoughPeaks_",step)

    return(list(DAT, dat))
}

checkLength <- function(...){

  nCompoundsOld <- data.table::uniqueN(processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE], by = c("ID"))
  nReplicatesOld <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("Batch"))
  nSeriesOld <- data.table::uniqueN(processingGroup[get(paste0("enoughPeaks_", step-1)) %in% TRUE], by = c("groupIndices"))

  nCompoundsNew <- data.table::uniqueN(processingGroup[get(paste0("enoughPeaks_", step)) %in% TRUE], by = c("ID"))
  nReplicatesNew <- data.table::uniqueN(processingFeature[color %in% "black"], by = c("Batch"))
  nSeriesNew <- data.table::uniqueN(processingGroup[get(paste0("enoughPeaks_", step)) %in% TRUE], by = c("groupIndices"))

  message(paste0("Removed ", nSeriesOld - nSeriesNew, " ", Series, " with less than ", MIN_FEATURE, " Signals.","\n",
                 "Remaining Data set with ", nCompoundsNew, " ", Compounds," and ", nReplicatesNew, " Batch(es) -> ", nSeriesNew, " ", Series),"\n--------------------------------------------------------\n")

  stopifnot(exprs = {
    "all Compounds were removed" = nCompoundsNew - data.table::uniqueN(processingGroup[get(paste0("enoughPeaks_", step)) %in% FALSE], by = c("ID")) > 0
    "all Dilution/Concentration-Series were removed" = nSeriesNew - data.table::uniqueN(processingGroup[get(paste0("enoughPeaks_", step)) %in% FALSE], by = c("groupIndices")) > 0
  })

  cutoff <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step)) %in% FALSE, groupIndices]]
  processingFeature <- processingFeature[groupIndices %in% processingGroup[get(paste0("enoughPeaks_", step)) %in% TRUE, groupIndices]]

  return(list(processingFeature, cutoff))


}

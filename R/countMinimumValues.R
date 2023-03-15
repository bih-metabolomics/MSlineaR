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
countMinimumValue <- function(DAT, MIN_FEATURE = 3, step){

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

  colnames(dat)[which(colnames(dat) %in% "N")] <- paste0("N_",step)
  colnames(dat)[which(colnames(dat) %in% "enoughPeaks")] <- paste0("enoughPeaks_",step)

    return(list(DAT, dat))
}

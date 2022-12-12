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
countMinimumValue <- function(DAT, MIN_FEATURE = 3, ...){

  stopifnot(exprs ={
    "Parameter MIN_FEATURE needs to be greater or equal than 3" = MIN_FEATURE >= 3


  })

  dat <- data.table::copy(DAT)

  dat <- unique(dat[, N := sum(!is.na(Y)),by = .(groupIndices)][ , ":=" ( enoughPeaks = ifelse( N >= MIN_FEATURE, TRUE, FALSE),
                                                     Comment = ifelse( N >= MIN_FEATURE, "EnoughPeaks", "notEnoughPeaks"))][
                                                       ,list(groupIndices,ID, REPLICATE,  N, enoughPeaks, Comment)])

    return(dat)
}

#' Remove signals below a blank-derived noise threshold
#'
#' @description
#' Filters signals based on their signal-to-blank ratio by comparing feature
#' intensities against the median intensity observed in blank samples.
#'
#' Signals with intensities less than or equal to a user-defined multiple of
#' the median blank intensity are considered indistinguishable from analytical
#' background noise and are excluded from further linearity assessment.
#'
#' @details
#' For each feature-specific dilution or concentration series, the median
#' intensity of the corresponding blank samples is calculated.
#'
#' Signals are flagged when:
#'
#' \deqn{
#' Signal \le MedianBlank \times noise
#' }
#'
#' Flagged observations are:
#' \itemize{
#'   \item Annotated in the `Comment` column
#'   \item Assigned the colour `"grey"`
#'   \item Excluded from subsequent analysis by setting `Y_sb` to `NA`
#' }
#'
#' If no blank intensity is available, all observations are retained and
#' annotated with `"noBlank"`.
#'
#' An additional continuity check is performed to remove isolated signal
#' segments separated from the main dilution-response trend after blank
#' filtering.
#'
#' @param dats A data.table containing a single feature-specific dilution or
#' concentration series.
#'
#' @param blanks A data.table containing blank measurements corresponding to
#' the same feature.
#'
#' @param y Character. Column name of the untransformed intensity variable.
#'
#' @param y_trans Character. Column name of the transformed intensity variable
#' used for downstream modelling.
#'
#' @param noise Numeric. Signal-to-blank threshold factor. Signals must exceed
#' `noise × median(blank)` to be retained.
#'
#' @return A data.table containing:
#' \describe{
#'   \item{signalBlankRatio}{Logical indicator of whether the signal failed
#'   the signal-to-blank criterion.}
#'
#'   \item{medBlank}{Median blank intensity used for filtering.}
#'
#'   \item{Y_sb}{Filtered response variable used for subsequent analysis.}
#' }
#'
#' @references
#' Signal-to-blank filtering is commonly applied in untargeted metabolomics
#' to remove analytical background signals and improve data quality before
#' statistical analysis.
#'
#' @keywords internal

#' @import data.table
#' @import dplyr
#' @importFrom stats median
trimm_signalBlank <- function(dats, blanks, y, y_trans, noise){

  data.table::setDT(dats)
  data.table::setDT(blanks)

  dat <- data.table::copy(dats)
  blank <- data.table::copy(blanks)
  blank <- dplyr::filter(blank, groupIndices %in% dat$groupIndices)
  dat$'signalBlankRatio' <- FALSE
  dat$Y_sb <- dat[[y_trans]]


  medblank <- median(blank[[y]], na.rm = T)

  if(!is.na(medblank)){




  dat$'signalBlankRatio'[dat[[y]] <= (medblank * noise)] <- TRUE
  dat$medBlank <- medblank
  dat$color[dat$'signalBlankRatio' %in% TRUE] <- "grey"
  dat$Comment[dat$'signalBlankRatio' %in% FALSE] <- paste0(dat$Comment[dat$'signalBlankRatio' %in% FALSE], "_>s/b")
  dat$Comment[dat$'signalBlankRatio' %in% TRUE] <- paste0(dat$Comment[dat$'signalBlankRatio' %in% TRUE], "_<s/b")
  dat$'signalBlankRatio'[is.na(dat[[y_trans]])] <- NA

  dat$Y_sb[dat$'signalBlankRatio' %in% TRUE] <- NA

  lost <- diff(dat[dat$color %in% "black", DilutionPoint]) > 2
  if(any(lost)){
    range <- dat[dat$color %in% "black", DilutionPoint]
    lost.pos <- which(lost)
    left <- length(1:lost.pos)
    right <- length((lost.pos + 1) : length(range))

    if(left < right & left < 3 ) {

      dat$'signalBlankRatio'[dat$color %in% "black"][1:lost.pos] <- TRUE


    } else if(right < left & right < 3 ){

      dat$'signalBlankRatio'[dat$color %in% "black"][ (lost.pos + 1) : length(lost)] <- TRUE

    }


    dat$color[dat$'signalBlankRatio' %in% TRUE] <- "grey"
    dat$Comment[dat$'signalBlankRatio' %in% TRUE] <- paste0(dat$Comment[dat$signalBlankRatio], "_<s/b")
    dat$Y_sb[dat$'signalBlankRatio' %in% TRUE] <- NA

  }

  }else{

    dat$Comment <- paste0(dat$Comment, "_noBlank")
    dat$medBlank <- NA

  }


  return(dat)

}







#' Remove non-linear boundaries and plateau regions from dilution-response curves
#'
#' @description
#' Identifies and removes non-linear regions at the lower and upper boundaries
#' of dilution or concentration series.
#'
#' The function iteratively trims observations until the remaining curve
#' exhibits the expected monotonic behaviour, where the lowest concentration
#' corresponds to the lowest response and the highest concentration corresponds
#' to the highest response.
#'
#' Additionally, plateau regions indicative of detector saturation or signal
#' compression are identified and removed using local slope analysis.
#'
#' @details
#' Non-linear behaviour frequently occurs at the extremes of the analytical
#' dynamic range due to:
#'
#' \itemize{
#'   \item Detector saturation at high concentrations
#'   \item Background noise at low concentrations
#'   \item Ion suppression effects
#'   \item Peak integration artefacts
#' }
#'
#' The trimming procedure consists of:
#'
#' \enumerate{
#'   \item Identification of upper non-linear regions.
#'   \item Identification of lower non-linear regions.
#'   \item Detection of plateau regions using local slope analysis.
#'   \item Annotation and exclusion of affected observations.
#' }
#'
#' Plateau detection is performed using the helper function
#' \code{findPlateaus()}.
#'
#' Observations identified as non-linear are excluded from downstream
#' modelling by setting `Y_trim` to `NA`.
#'
#' @param dats A data.table containing a single feature-specific dilution or
#' concentration series.
#'
#' @param y Character. Column name of the response variable used for trimming.
#'
#' @param x Character. Column name of the concentration or dilution variable.
#'
#' @param SLOPE_RATIO Numeric. Threshold used for plateau detection.
#' Slopes smaller than `SLOPE_RATIO × middleSlope` are considered part of a
#' plateau region. Default is `0.5`.
#'
#' @return A data.table containing:
#' \describe{
#'   \item{trim}{Logical indicator identifying trimmed observations.}
#'
#'   \item{flat_slope}{Logical indicator identifying regions with reduced slope.}
#'
#'   \item{Plateau}{Logical indicator identifying plateau regions.}
#'
#'   \item{Y_trim}{Response variable after trimming. Removed observations are
#'   set to `NA`.}
#' }
#'
#' @section Rationale:
#' Removing non-linear boundaries allows identification of the largest
#' feature-specific interval exhibiting approximately proportional
#' dilution-response behaviour, which is subsequently used for linearity
#' assessment.
#'
#' @keywords internal
#'
#' @import data.table
#' @import dplyr
#' @importFrom tidyr unite
trimEnds <- function(dats, y = parent.frame()$Y, x = parent.frame()$X, SLOPE_RATIO = 0.5){ # thresh=0

  dats$trim <- FALSE
  dats$flat_slope <- FALSE
  dats$Plateau <- FALSE
  dats$trim[is.na(dats[[y]])] <- NA # | dats$OutlierFOD %in% TRUE
  dat <- data.table::setorderv(dats,x)[!is.na(get(y))]# & !OutlierFOD %in% TRUE]

  dat$flat_slope <- findPlateaus(dats = dat,y = y, x = x, slope_ratio = SLOPE_RATIO)

  #browser()
  if (data.table::last(dat[[y]]) != max(dat[[y]])) {

    fin <- FALSE
    i = 1
    max_red <- list()
    while(fin == FALSE){
    maxPoint <- which.max(dat[[y]])
    maxminPoint <- which(dat[[y]] %in% min(dat[[y]][maxPoint:length(dat[[y]])]))
    dat.reduced.max <- data.table::copy(dat)[
      get(y) >= dat[[y]][maxminPoint],# + thresh,
      ':=' (trim = TRUE,
            Comment = "trim: >lastPoint")]
    # dat.reduced.max <- data.table::copy(dat)[
    #     which.max(get(tidyselect::all_of(y))) : data.table::last(get(tidyselect::all_of(y))),#get(tidyselect::all_of(y)) >= dat[,data.table::last(get(tidyselect::all_of(y)))] + thresh,
    #     ':=' (trim = TRUE,
    #           Comment = "trim: >lastPoint")]
    dat.reduced.max[nrow(dat.reduced.max),
                    ':=' (Comment = "trim: lastPoint",
                          trim = TRUE)]
if(any(dat.reduced.max$trim %in% FALSE)){
    if(data.table::last(dat.reduced.max[[y]][dat.reduced.max$trim %in% FALSE]) != max(dat.reduced.max[[y]][dat.reduced.max$trim %in% FALSE])){
      fin <- FALSE
      max_red[[i]] <- dat.reduced.max[dat.reduced.max$trim %in% TRUE,]
      dat <- dat.reduced.max[dat.reduced.max$trim %in% FALSE,]
      i = i +1

    } else {
      fin <- TRUE
      max_red[[i]] <- dat.reduced.max
      dat.reduced.max <- data.table::rbindlist(max_red) |> dplyr::arrange(DilutionPoint)

    }}else{fin <- TRUE
    max_red[[i]] <- dat.reduced.max
    dat.reduced.max <- data.table::rbindlist(max_red) |> dplyr::arrange(DilutionPoint)
        }

    }


  } else {

    dat.reduced.max <- dat[, trim := FALSE]

  }

  dat.reduced.max <- dat.reduced.max[flat_slope & rev(cumprod(rev(flat_slope))),
                         ':=' (trim = TRUE,
                               Plateau = TRUE,
                               Comment = paste0(Comment, "_trim: high Plateau"))]



  if (dat[[y]][1] != min(dat[[y]])){

    fin <- FALSE
    i = 1
    min_red <- list()
    while(fin == FALSE){

    minPoint <- which.min(dat[[y]])
    minmaxPoint <- which.max(dat[[y]][1:minPoint])
    dat.reduced.min <- data.table::copy(dat)[get(y) <= dat[,get(y)][minmaxPoint], #- thresh,
                                             ':=' (trim = TRUE,
                                                   Comment = "trim: <firstPoint")]
    # dat.reduced.min <- data.table::copy(dat)[1: which.min(get(tidyselect::all_of(y))), #get(tidyselect::all_of(y)) <= dat[,get(tidyselect::all_of(y))][1] - thresh,
    #                         ':=' (trim = TRUE,
    #                               Comment = "trim: <firstPoint")]
    dat.reduced.min[1, ':=' (Comment = "trim: firstPoint",
                             trim = TRUE)]
if(any(dat.reduced.min$trim %in% FALSE)){
    if(dat.reduced.min[[y]][dat.reduced.min$trim %in% FALSE][1] != min(dat.reduced.min[[y]][dat.reduced.min$trim %in% FALSE])){
      fin <- FALSE
      min_red[[i]] <- dat.reduced.min[dat.reduced.min$trim %in% TRUE,]
      dat <- dat.reduced.min[dat.reduced.min$trim %in% FALSE,]
      i = i +1

    } else {
      fin <- TRUE
      min_red[[i]] <- dat.reduced.min
      dat.reduced.min <- data.table::rbindlist(min_red) |> dplyr::arrange(DilutionPoint)

    }} else{fin <- TRUE}

    }



  } else {dat.reduced.min <- dat[, trim := FALSE]}

  dat.reduced.min <- dat.reduced.min[flat_slope & cumprod(flat_slope),
                                     ':=' (trim = TRUE,
                                           Plateau = TRUE,
                                           Comment = paste0(Comment, "_trim: low Plateau"))]

  tmp <-  dplyr::full_join(x = dat.reduced.min[trim %in% TRUE],
                           y = dat.reduced.max[trim %in% TRUE],
                           by = colnames(dat.reduced.max)[!(colnames(dat.reduced.max) %in% c("Comment", "flat_slope", "Plateau"))],
                           suffix = c(".x", ".y")) |>
    tidyr::unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE, sep = " / ") |>
    dplyr::mutate(Plateau = ifelse(as.logical(Plateau.x) + as.logical(Plateau.y) >=1,TRUE, FALSE)) |>
    dplyr::mutate(flat_slope =  ifelse(as.logical(flat_slope.x) + as.logical(flat_slope.y) >=1, TRUE, FALSE))

  data.table::setorder(tmp,DilutionPoint)
  data.table::setorder(dats,DilutionPoint)

  dats[IDintern %in% tmp$IDintern , Comment := paste(Comment,tmp$Comment, sep = "_")]
  dats[IDintern %in% tmp$IDintern, trim := ifelse(is.null(tmp$trim), FALSE, tmp$trim)]
  dats[IDintern %in% tmp$IDintern, flat_slope := ifelse(is.null(tmp$flat_slope), FALSE, tmp$flat_slope)]
  dats[IDintern %in% tmp$IDintern, Plateau := ifelse(is.null(tmp$Plateau), FALSE, tmp$Plateau)]
  dats$Comment[dats$trim %in% FALSE] <- paste0(dats$Comment[dats$trim %in% FALSE], "_NoTrim")
  dats$color[dats$trim %in% TRUE] <- "grey"
  dats$trim[is.na(dats[[y]])] <- NA
  dats$Y_trim <- dats[[y]]
  dats$Y_trim[dats$trim %in% TRUE] <- NA


  return(dats)
}



#' Detect plateau regions using local slope analysis
#'
#' @description
#' Identifies regions at the lower or upper end of a dilution-response curve
#' where the response changes more slowly than expected relative to the central
#' portion of the curve.
#'
#' Such regions are characteristic of detector saturation, signal compression,
#' or background noise and may indicate that observations lie outside the
#' useful analytical range.
#'
#' @details
#' The function:
#'
#' \enumerate{
#'   \item Determines the observation closest to half-maximal response.
#'   \item Estimates a local slope around this midpoint.
#'   \item Calculates slopes between all neighbouring observations.
#'   \item Compares local slopes against the midpoint slope.
#'   \item Flags observations where slopes fall below the specified fraction
#'   of the midpoint slope.
#' }
#'
#' An observation is classified as belonging to a plateau only if both
#' neighbouring slopes satisfy the plateau criterion.
#'
#' @param dats A data.table containing a single feature-specific dilution or
#' concentration series.
#'
#' @param y Character. Column name of the response variable.
#'
#' @param x Character. Column name of the concentration or dilution variable.
#'
#' @param slope_ratio Numeric. Fraction of the midpoint slope used as the
#' plateau threshold. Default is `0.5`.
#'
#' @return A logical vector with one value per observation indicating whether
#' the observation belongs to a plateau region.
#'
#' @section Interpretation:
#' A value of `TRUE` indicates that the local slope surrounding the
#' observation is substantially reduced relative to the central slope of the
#' curve and therefore may lie outside the useful analytical range.
#'
#' @keywords internal
findPlateaus <- function(dats, y, x , slope_ratio){

  dat <- dats
  # calculate point nearest half max intensity
  int50 <- DescTools::Closest(x = dat[[y]] ,a = (min(dat[[y]]) + max(dat[[y]]))/2, which = TRUE, na.rm = T)
  if(length(int50) > 1) int50 <- max(int50)
  if(int50 == length(dat[[x]])) int50 <- length(dat[[x]]) -1
  if(int50 == 1) int50 = 2

  #create linear regression line going through int50

  middleSlope <- suppressWarnings(summary(lm(dat[[y]][(int50 - 1) : (int50 + 1)] ~ dat[[x]][(int50 - 1) : (int50 + 1)]))$coefficient[2])
  allSlopes <- diff(dat[[y]])/diff(dat[[x]])

  allSlopes_a <- c(allSlopes, dplyr::last(allSlopes))
  allSlopes_b <- c(dplyr::first(allSlopes), allSlopes)


  plateau_a <- round(allSlopes_a,2) <= round(slope_ratio * middleSlope,2)
  plateau_b <- round(allSlopes_b,2) <= round(slope_ratio * middleSlope,2)

  plateau <- (plateau_a + plateau_b) == 2

  #if(plateau[1] %in% TRUE) {plateau <-  c(TRUE, plateau)} else {plateau <-  c(FALSE, plateau)}
  #plateau <- data.frame(DilutionPoint = dat$DilutionPoint, plateau = plateau)

  return(plateau)


  }

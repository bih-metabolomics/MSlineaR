#' Title
#'
#' @param dats
#' @param blanks
#' @param y
#' @param noise
#'
#' @return
#' @export
#'
#' @examples
trimm_signalBlank <- function(dats, blanks, y, y_trans,yblank, noise){

  setDT(dats)
  setDT(blanks)

  dat <- data.table::copy(dats)
  blank <- data.table::copy(blanks)
  blank <- dplyr::filter(blank, groupIndices %in% dat$groupIndices)
  dat$'signalBlankRatio' <- FALSE
  dat$Y_sb <- dat[[y_trans]]


  medblank <- median(blank[[yblank]], na.rm = T)

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







#' Title
#'
#' @param dats
#' @param y
#' @param x
#' @param thresh
#'
#' @return
#' @export
#'
#' @examples
trimEnds <- function(dats, y = parent.frame()$Y, x = parent.frame()$X, thresh=0){

  dats$trim <- FALSE
  dats$trim[is.na(dats[[y]])] <- NA # | dats$OutlierFOD %in% TRUE
  dat <- data.table::setorderv(dats,x)[!is.na(get(y))]# & !OutlierFOD %in% TRUE]

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
      dat.reduced.max <- rbindlist(max_red) |> dplyr::arrange(DilutionPoint)

    }}else{fin <- TRUE
    max_red[[i]] <- dat.reduced.max
    dat.reduced.max <- rbindlist(max_red) |> dplyr::arrange(DilutionPoint)
        }

    }


  } else {dat.reduced.max <- dat[, trim := FALSE]}

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
    if(data.table::last(dat.reduced.min[[y]][dat.reduced.min$trim %in% FALSE]) != min(dat.reduced.min[[y]][dat.reduced.min$trim %in% FALSE])){
      fin <- FALSE
      min_red[[i]] <- dat.reduced.min[dat.reduced.min$trim %in% TRUE,]
      dat <- dat.reduced.min[dat.reduced.min$trim %in% FALSE,]
      i = i +1

    } else {
      fin <- TRUE
      min_red[[i]] <- dat.reduced.min
      dat.reduced.min <- rbindlist(min_red) |> dplyr::arrange(DilutionPoint)

    }} else{fin <- TRUE}

    }



  } else {dat.reduced.min <- dat[, trim := FALSE]}

  tmp <-  dplyr::full_join(x = dat.reduced.min[trim %in% TRUE],
                           y = dat.reduced.max[trim %in% TRUE],
                           by = colnames(dat.reduced.max)[colnames(dat.reduced.max) != "Comment"],
                           suffix = c(".x", ".y")) |>
    tidyr::unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE, sep = " / ")
  data.table::setorder(tmp,DilutionPoint)
  data.table::setorder(dats,DilutionPoint)

  dats[IDintern %in% tmp$IDintern , Comment := paste(Comment,tmp$Comment, sep = "_")]
  dats[IDintern %in% tmp$IDintern, trim := ifelse(is.null(tmp$trim), FALSE, tmp$trim)]
  dats$Comment[dats$trim %in% FALSE] <- paste0(dats$Comment[dats$trim %in% FALSE], "_NoTrim")
  dats$color[dats$trim %in% TRUE] <- "grey"
  dats$trim[is.na(dats[[y]])] <- NA
  dats$Y_trim <- dats[[y]]
  dats$Y_trim[dats$trim %in% TRUE] <- NA


  return(dats)
}



#' Title
#'
#' @param dats
#' @param y
#' @param x
#'
#' @return
#' @export
#'
#' @examples
trim_pos_associated <- function(dats, y, x, MIN_Feature){
data.table::setDT(dats)

  dats$trimPos <- FALSE
  dats$trimPos[is.na(dats[[y]])] <- NA # | dats$OutlierFOD %in% TRUE
  dat <- data.table::setorderv(dats,x)[!is.na(get(y))]# & !OutlierFOD %in% TRUE]


  IsPositivAssociated =rle(c(dat[[y]][1] < dat[[y]][2], (dat[[y]][-1] - data.table::shift(dat[[y]], 1, type = "lag")[-1]) > 0))
  if(any(IsPositivAssociated$values %in% TRUE & IsPositivAssociated$lengths >= MIN_Feature)){

    if(length(max(IsPositivAssociated$lengths[IsPositivAssociated$values %in% TRUE])) == 1 ){
      PA_TRUE_length_max <- max(IsPositivAssociated$lengths[IsPositivAssociated$values %in% TRUE])
      PA_TRUE_length_max_pos <- which(IsPositivAssociated$lengths %in% PA_TRUE_length_max)

      rangeStart = ifelse(PA_TRUE_length_max_pos == 1, 1, cumsum(IsPositivAssociated$lengths)[PA_TRUE_length_max_pos-1] +1)

      PA_TRUE_range <- c(rangeStart:cumsum(IsPositivAssociated$lengths)[PA_TRUE_length_max_pos])

      exspected_min <- min(PA_TRUE_range)
      exspected_max <- max(PA_TRUE_range)

      if(dat[[y]][exspected_max] != data.table::last(dat[[y]]) &
         any(dat[[y]][exspected_max : which.max(dat[[y]])] < dat[[y]][exspected_max])){

        last_min <- which(dat[[y]] %in% min(dat[[y]][exspected_max:which.max(dat[[y]])]))
        dat <- dat[get(y) >= dat[last_min, get(y)],
                   ':=' (trimPos = TRUE,
                         Comment = paste(Comment,"trimPos", sep = "_"))]

      }else if(dat[[y]][exspected_min] != dat[[y]][1] &
               any(dat[[y]][1: exspected_min] > dat[[y]][exspected_min])){

        first_max <- which(dat[[y]] %in% max(dat[[y]][1: exspected_min]))

        second_max <- which(dat[[y]] %in% min(dat[[y]][dat[[y]] > dat[[y]][first_max]])) -1
        dat <- dat[1:second_max,
                        ':=' (trimPos = TRUE,
                              Comment = paste(Comment,"trimPos", sep = "_"))]
      }else{

        dat$trimPos[PA_TRUE_range] <- FALSE
        dat$trimPos[-PA_TRUE_range] <- TRUE




      }

      #if(PA_TRUE_length_max == )



    }else{

      dat$Comment = paste0(dat$Comment[dat$trimPos %in% FALSE], "_Multiple.PA.ranges")
    }


  }

  dats <- dplyr::full_join(dat, dats[!IDintern %in% dat$IDintern ], by = colnames(dats))

  dats$Comment[dats$trimPos %in% FALSE] <- paste0(dats$Comment[dats$trimPos %in% FALSE], "_NoTrimPos")
  dats$color[dats$trimPos %in% TRUE] <- "grey"
  dats$Y_trimPos <- dats[[y]]
  dats$Y_trimPos[dats$trimPos %in% TRUE] <- NA

  return(dats)


}





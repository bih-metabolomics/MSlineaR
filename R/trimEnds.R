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
trimEnds <- function(dats, y="YLog", x="XLog", thresh=0){

  dats$trim <- FALSE
  dats$trim[is.na(dats[, get(y)]) | dats$OutlierFOD %in% TRUE] <- NA
  dat <- data.table::setorder(dats,DilutionPoint)[!is.na(get(y)) & !OutlierFOD %in% TRUE]

  #browser()
  if (data.table::last(dat[[tidyselect::all_of(y)]]) != max(dat[[tidyselect::all_of(y)]])) {

    maxPoint <- which.max(dat[[tidyselect::all_of(y)]])
    maxminPoint <- which(dat[[tidyselect::all_of(y)]] %in% min(dat[[tidyselect::all_of(y)]][maxPoint:length(dat[[tidyselect::all_of(y)]])]))
    dat.reduced.max <- data.table::copy(dat)[
        get(tidyselect::all_of(y)) >= dat[,get(tidyselect::all_of(y))][maxminPoint],# + thresh,
        ':=' (trim = TRUE,
              Comment = "trim: >lastPoint")]
    # dat.reduced.max <- data.table::copy(dat)[
    #     which.max(get(tidyselect::all_of(y))) : data.table::last(get(tidyselect::all_of(y))),#get(tidyselect::all_of(y)) >= dat[,data.table::last(get(tidyselect::all_of(y)))] + thresh,
    #     ':=' (trim = TRUE,
    #           Comment = "trim: >lastPoint")]
    dat.reduced.max[nrow(dat.reduced.max),
                    ':=' (Comment = "trim: lastPoint",
                          trim = TRUE)]

  } else {dat.reduced.max <- dat[, trim := FALSE]}

  if (dat[[tidyselect::all_of(y)]][1] != min(dat[[tidyselect::all_of(y)]])){

    minPoint <- which.min(dat[[tidyselect::all_of(y)]])
    minmaxPoint <- which.max(dat[[tidyselect::all_of(y)]][1:minPoint])
    dat.reduced.min <- data.table::copy(dat)[get(tidyselect::all_of(y)) <= dat[,get(tidyselect::all_of(y))][minmaxPoint], #- thresh,
                             ':=' (trim = TRUE,
                                   Comment = "trim: <firstPoint")]
    # dat.reduced.min <- data.table::copy(dat)[1: which.min(get(tidyselect::all_of(y))), #get(tidyselect::all_of(y)) <= dat[,get(tidyselect::all_of(y))][1] - thresh,
    #                         ':=' (trim = TRUE,
    #                               Comment = "trim: <firstPoint")]
     dat.reduced.min[1, ':=' (Comment = "trim: firstPoint",
                             trim = TRUE)]



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

  return(dats)
}

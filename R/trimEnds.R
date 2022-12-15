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
  dats$Comment <- ""
  dats$trim[is.na(dats[, get(y)]) | dats$OutlierFOD %in% TRUE] <- NA
  dat <- setorder(dats,DilutionPoint)[!is.na(get(y)) & !OutlierFOD %in% TRUE]

  #browser()
  if (last(dat[[tidyselect::all_of(y)]]) != max(dat[[tidyselect::all_of(y)]])) {
    dat.reduced.max <- copy(dat)[
        get(tidyselect::all_of(y)) >= last(get(tidyselect::all_of(y))) - thresh,
        ':=' (trim = TRUE,
              Comment = "trim: >lastPoint")]
    dat.reduced.max[nrow(dat.reduced.max),
                    Comment := stringr::str_replace(Comment, "trim: >lastPoint", "trim: lastPoint")]

  } else {dat.reduced.max <- dat[, trim := FALSE]}

  if (dat[[tidyselect::all_of(y)]][1] != min(dat[[tidyselect::all_of(y)]])){

    dat.reduced.min <- copy(dat)[get(tidyselect::all_of(y)) <= get(tidyselect::all_of(y))[1] + thresh,
                            ':=' (trim = TRUE,
                                  Comment = "trim: <firstPoint")]
    dat.reduced.min[1, Comment := stringr::str_replace(Comment, "trim: <firstPoint", "trim: firstPoint")]

    } else {dat.reduced.min <- dat[, trim := FALSE]}

  tmp <-  dplyr::full_join(x = dat.reduced.min[trim %in% TRUE],
                    y = dat.reduced.max[trim %in% TRUE],
                    by = colnames(dat.reduced.max)[colnames(dat.reduced.max) != "Comment"],
                    suffix = c(".x", ".y")) |>
    tidyr::unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE, sep = " / ")
  setorder(tmp,DilutionPoint)
  setorder(dats,DilutionPoint)
  dats[IDintern %in% tmp$IDintern & (!Comment %in% c(NA, NULL, "", " ")), Comment := paste(Comment,tmp$Comment, sep = "_")]
  dats[IDintern %in% tmp$IDintern & (Comment %in% c(NA, NULL, "", " ")), Comment := tmp$Comment]
  dats[IDintern %in% tmp$IDintern, trim := ifelse(is.null(tmp$trim), FALSE, tmp$trim)]


  return(dats)
}

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
trimEnds <- function(dats, y="IntensityNorm", x="DilutionPoint", thresh=0){
  dats[ , Comment := as.character(Comment)]
  dat <- setorder(dats,DilutionPoint)[!is.na(get(y)) & outlierFOD %in% FALSE]
  dat[ , Comment := as.character(Comment)]


  #browser()
  if (last(dat[[all_of(y)]]) != max(dat[[all_of(y)]])) {
    dat.reduced.max <- copy(dat)[
      , Comment := as.character(Comment)][
        get(all_of(y)) >= last(get(all_of(y))) - thresh,
        ':=' (color = "grey", trim = TRUE,
              Comment = "trim: >lastPoint")]
    dat.reduced.max[nrow(dat.reduced.max),
                    Comment := str_replace(Comment, "trim: >lastPoint", "trim: lastPoint")]

  } else {dat.reduced.max <- dat[, trim := NA]}

  if (dat[[all_of(y)]][1] != min(dat[[all_of(y)]])){

    dat.reduced.min <- copy(dat)[
      , Comment := as.character(Comment)][get(all_of(y)) <= get(all_of(y))[1] + thresh,
                            ':=' (color = "grey", trim = TRUE,
                                  Comment = "trim: <firstPoint")]
    dat.reduced.min[1, Comment := str_replace(Comment, "trim: <firstPoint", "trim: firstPoint")]

    } else {dat.reduced.min <- dat[, trim := NA]}

  tmp <-  full_join(x = dat.reduced.min[trim %in% TRUE],
                    y = dat.reduced.max[trim %in% TRUE],
                    by = colnames(dat.reduced.max)[colnames(dat.reduced.max) != "Comment"],
                    suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE, sep = " / ")
  setorder(tmp,DilutionPoint)
  setorder(dats,DilutionPoint)
  dats[ , Comment := as.character(Comment)][IDintern %in% tmp$IDintern & (!Comment %in% c(NA, NULL, "", " ")), Comment := paste(Comment,tmp$Comment, sep = "_")]
  dats[ , Comment := as.character(Comment)][IDintern %in% tmp$IDintern & (Comment %in% c(NA, NULL, "", " ")), Comment := tmp$Comment]
  dats[IDintern %in% tmp$IDintern, ':=' (pch = tmp$pch, color = tmp$color, trim = ifelse(is.null(tmp$trim), NA, tmp$trim))]


  return(dats)
}

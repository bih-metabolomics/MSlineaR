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
  dat <- dats %>% drop_na(all_of(y)) %>% filter(outlierFOD %in% FALSE) %>% arrange(DilutionPoint) |>  dplyr::select(groupIndices, IDintern, all_of(x), all_of(y), Comment, color, pch)
  #browser()
  if (last(dat[[all_of(y)]]) != max(dat[[all_of(y)]])) {
    dat.reduced.max <- dat %>%
      filter(get(all_of(y)) >= last(get(all_of(y))) - thresh) %>%
      mutate(Comment = "trim:>lastPoint")
    dat.reduced.max$Comment[nrow(dat.reduced.max)] <- "trim:>lastPoint"
  } else {dat.reduced.max <- dat}

  if (dat[[all_of(y)]][1] != min(dat[[all_of(y)]])){
    dat.reduced.min <- dat %>%
      filter(get(all_of(y)) <= get(all_of(y))[1] + thresh) %>%
      mutate(Comment = "trim:<firstPoint")
    dat.reduced.min$Comment[1] <- "trim:firstPoint"
  } else {dat.reduced.min <- dat}

  tmp <-  full_join(x = dat.reduced.min,
                    y = dat.reduced.max,
                    by = colnames(dat.reduced.max)[colnames(dat.reduced.max) != "Comment"],
                    suffix = c(".x", ".y")) |>
    unite(Comment, c(Comment.x,Comment.y), remove = TRUE, na.rm = TRUE)

  tmp$color[str_detect(tmp$Comment, "trim:>lastPoint|trim:<firstPoint|trim:>lastPoint|trim:firstPoint")] <-  "grey"
  tmp$pch[str_detect(tmp$Comment, "trim:>lastPoint|trim:<firstPoint|trim:>lastPoint|trim:firstPoint")] <- 19
  if(any((!dats$IDintern %in% tmp$IDintern))){
    tmp <-  bind_rows( dats %>% filter(!IDintern %in% tmp$IDintern) |> dplyr::select(groupIndices, IDintern, Comment, color, pch),
                       tmp)}

  return(tmp |> dplyr::select(-all_of(x), -all_of(y)))
}

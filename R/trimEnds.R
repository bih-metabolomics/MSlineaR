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
trimEnds <- function(dats, y="IntensityNorm", x="DP", thresh=0){
  dat <- dats %>% drop_na(y) %>% filter(color != "red") %>% arrange(DilutionPoint)
  #browser()
  if (last(dat[[y]]) != max(dat[[y]])) {
    dat.reduced.max <- dat %>%
      filter(get(y) >= last(get(y)) - thresh) %>%
      mutate(Comment = "bigger")
    dat.reduced.max$Comment[nrow(dat.reduced.max)] <- "last"
  } else {dat.reduced.max <- dat}

  if (dat[[y]][1] != min(dat[[y]])){
    dat.reduced.min <- dat %>%
      filter(get(y) <= get(y)[1] + thresh) %>%
      mutate(Comment = "smaller")
    dat.reduced.min$Comment[1] <- "first"
  } else {dat.reduced.min <- dat}

  tmp <- full_join(dat.reduced.max,dat.reduced.min,
                   by = names(dat)[-which(arr.ind = T,names(dat) == "Comment" )]) %>%
    unite(Comment, c("Comment.x","Comment.y"), remove = TRUE, na.rm = TRUE)

  tmp$color[str_detect(tmp$Comment, "bigger|smaller|last|first")] <-  "grey"
  tmp$pch[str_detect(tmp$Comment, "bigger|smaller|last|first")] <- 19
  if(any((!dats$IDintern %in% tmp$IDintern))){
    tmp <-  full_join( dats %>% filter(!IDintern %in% tmp$IDintern), tmp, by = names(dat)[names(dats) %in% names(tmp)])}

  return(tmp)
}

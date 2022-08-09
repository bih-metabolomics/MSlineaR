#' Title
#'
#' @param dat
#' @param y
#' @param x
#' @param model
#'
#' @return
#' @export
#'
#' @examples
#'
chooseModel <- function(dat,
                        y="IntensityNorm",
                        x="DilutionPoint",
                        model=c("logistic", "linear", "quadratic")){

  dat <- dat %>% drop_na(all_of(y))
  dat <- dat %>% arrange(DilutionPoint)

  if ("logistic" %in% model) {
    logistic <- drc::drm(get(all_of(y)) ~ get(x), fct = L.3(), data = dat)
    cor.logistic <- cor(dat[all_of(y)], predict(logistic))
  } else{cor.logistic <- NA}

  if ("linear" %in% model) {
    linear <- lm(get(all_of(y)) ~ get(x), data = dat)
    cor.linear <- cor(dat[all_of(y)], predict(linear))
  } else{cor.linear <- NA}


  if ("quadratic" %in% model) {
    quadratic <- lm(get(all_of(y)) ~ poly(get(x), 2, raw = TRUE), data = dat)
    cor.quadratic <- cor(dat[all_of(y)], predict(quadratic))
  } else{
    cor.quadratic <- NA
  }


  # calculate correlation
  cor.poly <- c(
    "cor.logistic" = cor.logistic,
    "cor.linear" = cor.linear,
    "cor.quadratic" = cor.quadratic
  )

  cor.max <-  names(which(cor.poly == max(cor.poly, na.rm = T)))

  if ("cor.linear" %in% cor.max) {cor.max = "cor.linear"} else if (cor.max %in% c("cor.logistic", "cor.quadratic")) {cor.max = "cor.logistic"}  # if same correlation
  model <- get(substr(cor.max,5, 100 ))


  tmp <- list(
    "tmp" = list(
      "model.name" = substr(cor.max,5, 100 ),
      "model" = model,
      "cor" = max(cor.poly, na.rm = T)
    )
  )

  names(tmp) <- unique(dat$groupIndices)

  #SSE <- sum((fitted(cor.max) - dat$Intensity_norm)^2)

  return(tmp)

}

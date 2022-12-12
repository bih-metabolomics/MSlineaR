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
                        y="Intensity",
                        x="Concentration",
                        model=c("logistic", "linear", "quadratic")){

  dat <- setorder(dat,DilutionPoint)[!is.na(get(tidyselect::all_of(y)))]

  if ("logistic" %in% model) {
    logistic <- drc::drm(get(tidyselect::all_of(y)) ~ get(x), fct = drc::L.3(), data = dat)
    cor.logistic <- cor(dat[[tidyselect::all_of(y)]], predict(logistic))
  } else{cor.logistic <- NA}

  if ("linear" %in% model) {
    linear <- lm(get(tidyselect::all_of(y)) ~ get(x), data = dat)
    cor.linear <- cor(dat[[tidyselect::all_of(y)]], predict(linear))
  } else{cor.linear <- NA}


  if ("quadratic" %in% model) {
    quadratic <- lm(get(tidyselect::all_of(y)) ~ poly(get(x), 2, raw = TRUE), data = dat)
    cor.quadratic <- cor(dat[[tidyselect::all_of(y)]], predict(quadratic))
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
  #rm(list = c("logistic", "linear", "quadratic")[!c("logistic", "linear", "quadratic") %in% substr(cor.max,5, 100 )])


  if ("cor.linear" %in% cor.max) {cor.max = "cor.linear"} else if (all(c("cor.logistic", "cor.quadratic") %in% cor.max)) {cor.max = "cor.logistic"}  # if same correlation
  Model <- get(substr(cor.max,5, 100 ))

  Model = list("fit" = fitted(Model), "coefficients" = coef(Model), "residuals" = residuals(Model))



  tmp <- list(
    "tmp" = list(
      "model.name" = substr(cor.max,5, 100 ),
      "model" = Model,
      "cor" = max(cor.poly, na.rm = T)
    )
  )

  names(tmp) <- unique(dat$groupIndices)

  #SSE <- sum((fitted(cor.max) - dat$Intensity_norm)^2)
  return(tmp)

  rm(list = c("dat","cor.max", "Model", "cor.poly", "logistic", "linear" ))
  gc()

}

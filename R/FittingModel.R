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
                        model=c("logistic", "linear", "quadratic"),
                        SRES = 2){

  dat <- setorder(dat,DilutionPoint)[!is.na(get(tidyselect::all_of(y)))]
  dat$outlier <- FALSE

  if ("logistic" %in% model) {
    logistic <- drc::drm(get(tidyselect::all_of(y)) ~ get(x), fct = drc::L.3(), data = dat)
    if(any(residuals(logistic)/sd(residuals(logistic)) > SRES)){
      datOut <- dat
      datOut[[y]][which(residuals(logistic)/sd(residuals(logistic)) > SRES)] <- NA
      logisticOut <- drc::drm(get(tidyselect::all_of(y)) ~ get(x), fct = drc::L.3(), data = datOut)
      suppressWarnings(aic <- AIC(logistic, logisticOut))
      aic$AIC[is.infinite(aic$AIC)] <- -1000
      logistic1 <- get(row.names(aic)[which(aic$AIC %in% min(aic$AIC))])
    } else{
      logistic1 <- logistic
    }



    #cor.logistic <- cor(dat[[tidyselect::all_of(y)]], predict(logistic))
  } else{logistic1 <- NA}

  if ("linear" %in% model) {
    linear <- lm(get(tidyselect::all_of(y)) ~ get(x), data = dat)
    if(any(residuals(linear)/sd(residuals(linear)) > SRES)){
      datOut <- dat
      datOut[[y]][which(residuals(linear)/sd(residuals(linear)) > SRES)] <- NA
      linearOut <- lm(get(tidyselect::all_of(y)) ~ get(x), data = datOut)
      suppressWarnings(aic <- AIC(linear, linearOut))
      aic$AIC[is.infinite(aic$AIC)] <- -1000
      linear1 <- get(row.names(aic)[which(aic$AIC %in% min(aic$AIC))])
    } else{
      linear1 <- linear
    }

    #cor.linear <- cor(dat[[tidyselect::all_of(y)]], predict(linear))
  } else{linear1 <- NA}


  if ("quadratic" %in% model) {
    quadratic <- lm(get(tidyselect::all_of(y)) ~ poly(get(x), 2, raw = TRUE), data = dat)
    if(any(residuals(quadratic)/sd(residuals(quadratic)) > SRES)){
      datOut <- dat
      datOut[[y]][which(residuals(quadratic)/sd(residuals(quadratic)) > SRES)] <- NA
      quadraticOut <- lm(get(tidyselect::all_of(y)) ~ poly(get(x), 2, raw = TRUE), data = datOut)
      suppressWarnings(aic <- AIC(quadratic, quadraticOut))
      aic$AIC[is.infinite(aic$AIC)] <- -1000
      quadratic1 <- get(row.names(aic)[which(aic$AIC %in% min(aic$AIC))])
    } else{
      quadratic1 <- quadratic
    }



    #cor.quadratic <- cor(dat[[tidyselect::all_of(y)]], predict(quadratic))
  } else{
    quadratic1 <- NA
  }

  # select Model by using akaike information criterium

  suppressWarnings(aic <- AIC(linear1, logistic1, quadratic1))
  aic$AIC[is.infinite(aic$AIC)] <- -1000
  ModelName <- row.names(aic)[which(aic$AIC %in% min(aic$AIC))]
  if ("linear1" %in% ModelName) {ModelName = "linear1"} else if (all(c("logistic1", "quadratic1") %in% ModelName)) {ModelName = "logistic1"}  # if same correlation

  dat$outlier[which(residuals(get(gsub(pattern = "1", x = ModelName, replacement = "")))/sd(residuals(get(gsub(pattern = "1", x = ModelName, replacement = "")))) > SRES)] <- TRUE
  dat$outlier[as.integer(names(which(residuals(get(ModelName))/sd(residuals(get(ModelName))) > SRES)))] <- TRUE


  ## calculate correlation
  # cor.poly <- c(
  #   "cor.logistic" = cor.logistic,
  #   "cor.linear" = cor.linear,
  #   "cor.quadratic" = cor.quadratic
  # )
  # cor.max <-  names(which(cor.poly == max(cor.poly, na.rm = T)))
  #rm(list = c("logistic", "linear", "quadratic")[!c("logistic", "linear", "quadratic") %in% substr(cor.max,5, 100 )])



  Model <- ModelName

  Model = list("fit" = fitted(get(Model)), "coefficients" = coef(get(Model)), "residuals" = residuals(get(Model))/sd(residuals(get(Model))))



  tmp <- list(
    "tmp" = list(
      "model.name" = gsub(pattern = "1", x = ModelName, replacement = ""),
      "model" = Model,
      "aic" = aic$AIC[ row.names(aic) %in% ModelName],
      dat = dat
    )
  )

  names(tmp) <- unique(dat$groupIndices)

  #SSE <- sum((fitted(cor.max) - dat$Intensity_norm)^2)
  return(tmp)

  rm(list = c("dat","ModelName", "Model", "aic", "logistic", "linear" ))
  gc()

}

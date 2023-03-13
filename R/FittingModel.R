#' Title
#'
#' @param dats
#' @param y
#' @param x
#' @param model
#'
#' @return
#' @export
#' @import data.table
#' @examples
#'
chooseModel <- function(dats,
                        y = "Intensity",
                        x = "Concentration",
                        model = c("logistic", "linear", "quadratic"),
                        SDRES_MIN = 1,
                        SDRES = 2,
                        abbr,
                        R2_MIN = 0.95,
                        ...){
 # .datatable.aware=TRUE
  data.table::setDT(dats)
  dat <- data.table::copy(dats)[dats$color %in% "black",]
  outlierName <- paste0("Outlier", abbr)
  dat <- data.table::setorder(dat,DilutionPoint)[!is.na(get(y))]
  dat[[outlierName]] <- FALSE

  if ("logistic" %in% model) {
    logistic <- drc::drm(get(y) ~ get(x), fct = drc::L.3(), data = dat)
    if(any(abs(residuals(logistic)/sd(residuals(logistic))) > SDRES) & abs(sd(residuals(logistic))) > SDRES_MIN){
      datOutLog <- dat
      datOutLog[[y]][which(abs(residuals(logistic)/sd(residuals(logistic))) > SDRES)] <- NA
      logisticOut <- drc::drm(get(y) ~ get(x), fct = drc::L.3(), data = datOutLog)
      #RMSE
      logisticRMSE <- Metrics::rmse(dat[[y]], predict(logistic))
      logisticOutRMSE <- Metrics::rmse(na.exclude(datOutLog[[y]]), predict(logisticOut))
      logistic1 <- get(c("logistic", "logisticOut")[which(c(logisticRMSE, logisticOutRMSE) %in% min(logisticRMSE, logisticOutRMSE))])
    } else{
      logistic1 <- logistic
    }

    logistic1RMSE <- Metrics::rmse(logistic1$data$`get(y)`, predict(logistic1))

    #cor.logistic <- cor(dat[[tidyselect::all_of(y)]], predict(logistic))
  } else{
    logistic1 <- NA
    logistic1RMSE <- NA
  }

  if ("linear" %in% model) {
    linear <- lm(get(y) ~ get(x), data = dat)
    if(any(abs(residuals(linear)/sd(residuals(linear))) > SDRES) & abs(sd(residuals(linear))) > SDRES_MIN){
      datOutLin <- dat
      datOutLin[[y]][which(abs(residuals(linear)/sd(residuals(linear))) > SDRES)] <- NA
      linearOut <- lm(get(y) ~ get(x), data = datOutLin)
      #RMSE
      linearRMSE <- Metrics::rmse(dat[[y]], predict(linear))
      linearOutRMSE <- Metrics::rmse(na.exclude(datOutLin[[y]]), predict(linearOut))
      linear1 <- get(c("linear", "linearOut")[which(c(linearRMSE, linearOutRMSE) %in% min(linearRMSE, linearOutRMSE))])

    } else{
      linear1 <- linear
    }
    linear1RMSE <- Metrics::rmse(linear1$model$`get(y)`, predict(linear1))

    #cor.linear <- cor(dat[[tidyselect::all_of(y)]], predict(linear))
  } else{
    linear1 <- NA
    linear1RMSE <- NA
    }


  if ("quadratic" %in% model) {
    quadratic <- lm(get(y) ~ poly(get(x), 2, raw = TRUE), data = dat)
    if(any(abs(residuals(quadratic)/sd(residuals(quadratic))) > SDRES) & abs(sd(residuals(quadratic))) > SDRES_MIN){
      datOutQuad <- dat
      datOutQuad[[y]][which(abs(residuals(quadratic)/sd(residuals(quadratic))) > SDRES)] <- NA
      quadraticOut <- lm(get(y) ~ poly(get(x), 2, raw = TRUE), data = datOutQuad)
      #RMSE
      quadraticRMSE <- Metrics::rmse(dat[[y]], predict(quadratic))
      quadraticOutRMSE <- Metrics::rmse(na.exclude(datOutQuad[[y]]), predict(quadraticOut))
      quadratic1 <- get(c("quadratic", "quadraticOut")[which(c(quadraticRMSE, quadraticOutRMSE) %in% min(quadraticRMSE, quadraticOutRMSE))])
    } else{
      quadratic1 <- quadratic
    }

    quadratic1RMSE <- Metrics::rmse(quadratic1$model$`get(y)`, predict(quadratic1))

    #cor.quadratic <- cor(dat[[tidyselect::all_of(y)]], predict(quadratic))
  } else{
    quadratic1 <- NA
    quadratic1RMSE <- NA
  }

  # select Model by using RMSE
 # logistic1RMSE <- Metrics::rmse(logistic1$data$`get(tidyselect::all_of(y))`, predict(logistic1))
 # linear1RMSE <- Metrics::rmse(linear1$model$`get(tidyselect::all_of(y))`, predict(linear1))
 # quadratic1RMSE <- Metrics::rmse(quadratic1$model$`get(tidyselect::all_of(y))`, predict(quadratic1))

  ModelName <- c("logistic1", "linear1", "quadratic1")[which(round(c(logistic1RMSE, linear1RMSE, quadratic1RMSE),2) %in% min(round(c(logistic1RMSE, linear1RMSE, quadratic1RMSE),2),na.rm = TRUE))]
  if ("linear1" %in% ModelName) {ModelName = "linear1"} else if (all(c("logistic1", "quadratic1") %in% ModelName)) {ModelName = "logistic1"}  # if same correlation

 if(abs(sd(residuals(get(gsub(pattern = "1", x = ModelName, replacement = ""))))) > SDRES_MIN) dat[[outlierName]][which(abs(residuals(get(gsub(pattern = "1", x = ModelName, replacement = "")))/sd(residuals(get(gsub(pattern = "1", x = ModelName, replacement = ""))))) > SDRES)] <- TRUE
 if(abs(sd(residuals(get(ModelName)))) > SDRES_MIN & any(abs(sd(residuals(get(ModelName)))/sd(residuals(get(ModelName)))) > SDRES)){
   dat[IDintern %in% get(get(ModelName)$call$data)[!is.na(get(y))][which(abs(residuals(get(ModelName))/sd(residuals(get(ModelName)))) > SDRES)]$IDintern][[outlierName]] <- TRUE
}
  dat <- data.table::setorder(dplyr::full_join(dat, dats[!IDintern %in% dat$IDintern], by = colnames(dats)), DilutionPoint)
  dat$color[dat[[outlierName]] %in% TRUE] <- "red"
  dat$Comment[dat[[outlierName]] %in% TRUE] <- paste0(dat$Comment[dat[[outlierName]] %in% TRUE], "_Outlier",abbr)
  dat$Comment[dat[[outlierName]] %in% FALSE] <- paste0(dat$Comment[dat[[outlierName]] %in% FALSE], "_NoOutlier",abbr)



  modelNew <- if(ModelName %in% "logistic1"){drc::drm(get(y) ~ get(x), fct = drc::L.3(), data = dat[color %in% "black",])
  } else if(ModelName %in% "linear1"){lm(get(y) ~ get(x), data = dat[color %in% "black",])
      } else if(ModelName %in% "quadratic1"){lm(get(y) ~ poly(get(x), 2, raw = TRUE), data = dat[color %in% "black",])}

  ## calculate correlation
  # cor.poly <- c(
  #   "cor.logistic" = cor.logistic,
  #   "cor.linear" = cor.linear,
  #   "cor.quadratic" = cor.quadratic
  # )
  # cor.max <-  names(which(cor.poly == max(cor.poly, na.rm = T)))
  #rm(list = c("logistic", "linear", "quadratic")[!c("logistic", "linear", "quadratic") %in% substr(cor.max,5, 100 )])



  Model <- ModelName

  Model = list("fit" = fitted(get(Model)), "coefficients" = coef(get(Model)), "std.residuals" = residuals(get(Model))/sd(residuals(get(Model))), "RMSE" = c(c("logistic" = logistic1RMSE, "linear" = linear1RMSE, "quadratic" = quadratic1RMSE)))



  tmp <- list(
    "tmp" = list(
      "model.name" = gsub(pattern = "1", x = ModelName, replacement = ""),
      "model" = Model,
      "R2" = cor(dat[color %in% "black",][[y]], predict(modelNew))^2,
      "aboveMinCor" = cor(dat[color %in% "black",][[y]], predict(modelNew))^2 > R2_MIN,
      #"aic" = aic$AIC[ row.names(aic) %in% ModelName],
      dat = dat
    )
)
  names(tmp) <- unique(dat$groupIndices)

  #SSE <- sum((fitted(cor.max) - dat$Intensity_norm)^2)
  return(tmp)

  #rm(list = c("dat","ModelName", "Model", "aic", "logistic", "linear" ))
  #gc()

}

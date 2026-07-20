#' Detect and remove outliers using model-based regression fitting
#'
#' @description
#' `chooseModel()` identifies potential outliers within a dilution or
#' concentration series using a model-based approach.
#'
#' Three regression models can be evaluated:
#'
#' * Linear
#' * Quadratic
#' * Three-parameter logistic
#'
#' For each selected model:
#'
#' 1. The model is fitted to the complete series.
#' 2. Standardized residuals are calculated.
#' 3. Signals exceeding the specified residual threshold (`STDRES`) are
#'    considered potential outliers.
#' 4. A second model is fitted after excluding the potential outliers.
#' 5. The original and outlier-filtered models are compared.
#'
#' The final model is selected based on the lowest Bayesian Information
#' Criterion (BIC). Outliers identified by the selected model are flagged and
#' excluded from subsequent processing steps.
#'
#' To avoid excessive outlier removal in nearly perfect dilution series,
#' optional safeguards based on residual standard deviation (`SDRES_MIN`) and
#' local slope direction are applied.
#'
#' @param dats Data table containing a single dilution or concentration series
#'   after preprocessing and trimming.
#' @param y Character string specifying the column containing the transformed
#'   response variable (e.g. log-transformed signal intensity).
#' @param x Character string specifying the column containing the transformed
#'   dilution or concentration values.
#' @param model Character vector defining which regression models should be
#'   evaluated. Supported values are:
#'   \describe{
#'     \item{"linear"}{Simple linear regression}
#'     \item{"quadratic"}{Second-order polynomial regression}
#'     \item{"logistic"}{Three-parameter logistic regression}
#'   }
#' @param SDRES_MIN Numeric value specifying the minimum residual standard
#'   deviation required before residual-based outlier removal is performed.
#'   This prevents removal of points in nearly noise-free curves.
#' @param STDRES Numeric threshold for standardized residuals. Data points with
#'   absolute standardized residuals greater than this value are considered
#'   potential outliers.
#' @param abbr Character string used to label the outlier detection step
#'   (typically `"FOD"` for first outlier detection or `"SOD"` for second
#'   outlier detection).
#'
#' @return
#' A named list containing one entry per dilution/concentration series.
#' Each entry contains:
#'
#' \describe{
#'   \item{model.name}{Selected regression model.}
#'   \item{model}{Model summary information including fitted values,
#'   coefficients, standardized residuals and BIC values.}
#'   \item{dat}{Input data with additional columns indicating detected
#'   outliers and filtered response values.}
#' }
#'
#' Additional columns added to `dat` include:
#'
#' \describe{
#'   \item{OutlierFOD / OutlierSOD}{Logical indicator for detected outliers.}
#'   \item{Y_FOD / Y_SOD}{Response variable after removing detected outliers.}
#'   \item{Comment}{Annotation describing outlier status.}
#'   \item{color}{Visual status indicator used internally for plotting.}
#' }
#'
#' @details
#' The function evaluates whether removing potential outliers improves model
#' performance. Rather than relying solely on residual thresholds, the final
#' decision is based on model quality. This reduces the risk of incorrectly
#' removing biologically meaningful observations.
#'
#' When multiple models perform similarly, preference is given to linear
#' models, followed by logistic models.
#'
#' @references
#' Bayesian Information Criterion (BIC):
#' Schwarz G. (1978) Estimating the Dimension of a Model.
#' The Annals of Statistics 6(2):461-464.
#'
#' @seealso
#' \code{\link{trimEnds}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- chooseModel(
#'   dats = feature_data,
#'   y = "Y_trim",
#'   x = "X_trans",
#'   model = c("linear", "quadratic", "logistic"),
#'   STDRES = 2,
#'   SDRES_MIN = 0.1,
#'   abbr = "FOD"
#' )
#' }

chooseModel <- function(dats,
                        y = parent.frame()$Y,
                        x = parent.frame()$X,
                        model = c("logistic", "linear", "quadratic"),
                        SDRES_MIN = 0,
                        STDRES = 2,
                        abbr){
 # .datatable.aware=TRUE
  data.table::setDT(dats)
  dat <- data.table::copy(dats)[dats$color %in% "black",]
  outlierName <- paste0("Outlier", abbr)
  outlierY <- paste0("Y_", abbr)
  dat <- data.table::setorder(dat,DilutionPoint)[!is.na(get(y))]
  dat[[outlierName]] <- FALSE
  dat[[outlierY]] <- dat[[y]]

  if ("logistic" %in% model) {
    logistic <- drc::drm(get(outlierY) ~ get(x), fct = drc::L.3(), data = dat)
    logisticRMSE <- Metrics::rmse(dat[[outlierY]], predict(logistic))

    if(any(abs(residuals(logistic, typeRes = "standard")) > STDRES) ){ #& abs(sd(residuals(logistic))) > SDRES_MIN
      if(abs(sd(residuals(logistic))) > SDRES_MIN){
      datOutLog <- dat
      datOutLog[[outlierY]][which(abs(residuals(logistic, typeRes = "standard")) > STDRES)] <- NA
      datOutLog[[outlierName]][which(abs(residuals(logistic, typeRes = "standard")) > STDRES)] <- TRUE

      logisticOut <- drc::drm(get(outlierY) ~ get(x), fct = drc::L.3(), data = datOutLog)
      #RMSE
      logisticRMSE <- Metrics::rmse(dat[[outlierY]], predict(logistic))
      logisticOutRMSE <- Metrics::rmse(na.exclude(datOutLog[[outlierY]]), predict(logisticOut))
      logistic1 <- get(c("logistic", "logisticOut")[which(c(logisticRMSE, logisticOutRMSE) %in% min(logisticRMSE, logisticOutRMSE))])
      logistic1.dat <- datOutLog
      } else{
        refslopemin <- coef(logistic)[2]*100/3
        refslopemax <- coef(logistic)[2]*100*3
        slopes <- sapply(1:(nrow(dat)-1), function(i) coef(lm(dat = dat[i:(i+1)], get(outlierY) ~ get(x)))[2]*100)
        if(slopes[1] > refslopemin){ slopes <- c(refslopemin*2, slopes)} else{slopes <- c(slopes[1],refslopemin*2, slopes[-1])}

        pos_slope <- c(which(abs(residuals(logistic, typeRes = "standard")) > STDRES), which(abs(residuals(logistic, typeRes = "standard")) > STDRES) + 1)
        pos_slope <- pos_slope[pos_slope <= nrow(dat)]

        if(any(slopes[pos_slope] < 0)){

          datOutLog <- dat
          datOutLog[[outlierY]][which(abs(residuals(logistic, typeRes = "standard")) > STDRES)] <- NA
          datOutLog[[outlierName]][which(abs(residuals(logistic, typeRes = "standard")) > STDRES)] <- TRUE

          logisticOut <- drc::drm(get(outlierY) ~ get(x), fct = drc::L.3(), data = datOutLog)
          #RMSE
          logisticRMSE <- Metrics::rmse(dat[[outlierY]], predict(logistic))
          logisticOutRMSE <- Metrics::rmse(na.exclude(datOutLog[[outlierY]]), predict(logisticOut))
          logistic1 <- get(c("logistic", "logisticOut")[which(c(logisticRMSE, logisticOutRMSE) %in% min(logisticRMSE, logisticOutRMSE))])
          logistic1.dat <- datOutLog


        }else{

          logistic1 <- logistic
          logistic1.dat <- dat
        }

      }
    }else{
      logistic1 <- logistic
      logistic1.dat <- dat
    }

    logistic1RMSE <- logisticRMSE #Metrics::rmse(logistic1$data$`get(outlierY)`, predict(logistic1))

    #cor.logistic <- cor(dat[[tidyselect::all_of(y)]], predict(logistic))
  } else{
    logistic1 <- NA
    logistic1RMSE <- NA
  }

  if ("linear" %in% model) {
    linear <- lm(get(outlierY) ~ get(x), data = dat)
    linearRMSE <- Metrics::rmse(dat[[outlierY]], predict(linear))
    abs_std_residuals <- abs(rstandard(linear))

    if(any(abs_std_residuals > STDRES) ){ #& abs(sd(residuals(linear))) > SDRES_MIN
      if(abs(sd(residuals(linear))) > SDRES_MIN){
        datOutLin <- dat
        datOutLin[[outlierY]][which(abs_std_residuals > STDRES)] <- NA
        datOutLin[[outlierName]][which(abs_std_residuals > STDRES)] <- TRUE

        linearOut <- lm(get(outlierY) ~ get(x), data = datOutLin)
        #RMSE
        linearRMSE <- Metrics::rmse(dat[[outlierY]], predict(linear))
        linearOutRMSE <- Metrics::rmse(na.exclude(datOutLin[[outlierY]]), predict(linearOut))
        linear1 <- get(c("linear", "linearOut")[which(c(linearRMSE, linearOutRMSE) %in% min(linearRMSE, linearOutRMSE))])
        linear1.dat <- datOutLin
      } else{
        refslopemin <- coef(linear)[2]*100/3
        refslopemax <- coef(linear)[2]*100*3
        slopes <- sapply(1:(nrow(dat)-1), function(i) coef(lm(dat = dat[i:(i+1)], get(outlierY) ~ get(x)))[2]*100)
        if(slopes[1] > refslopemin){ slopes <- c(refslopemin*2, slopes)} else{slopes <- c(slopes[1],refslopemin*2, slopes[-1])}

        pos_slope <- c(which(abs_std_residuals > STDRES), which(abs_std_residuals > STDRES) + 1)
        pos_slope <- pos_slope[pos_slope <= nrow(dat)]

        if(any(slopes[pos_slope] < 0)){

          datOutLin <- dat
          datOutLin[[outlierY]][which(abs_std_residuals > STDRES)] <- NA
          datOutLin[[outlierName]][which(abs_std_residuals > STDRES)] <- TRUE

          linearOut <- lm(get(outlierY) ~ get(x), data = datOutLin)
          #RMSE
          linearRMSE <- Metrics::rmse(dat[[outlierY]], predict(linear))
          linearOutRMSE <- Metrics::rmse(na.exclude(datOutLin[[outlierY]]), predict(linearOut))
          linear1 <- get(c("linear", "linearOut")[which(c(linearRMSE, linearOutRMSE) %in% min(linearRMSE, linearOutRMSE))])
          linear1.dat <- datOutLin


        } else{

          linear1 <- linear
          linear1.dat <- dat
        }



      }
    }else{

      linear1 <- linear
      linear1.dat <- dat
    }
    linear1RMSE <- linearRMSE #Metrics::rmse(linear1$model$`get(outlierY)`, predict(linear1))

    #cor.linear <- cor(dat[[tidyselect::all_of(y)]], predict(linear))
  } else{
    linear1 <- NA
    linear1RMSE <- NA
    }


  if ("quadratic" %in% model) {
    quadratic <- lm(get(outlierY) ~ poly(get(x), 2, raw = TRUE), data = dat)
    quadraticRMSE <- Metrics::rmse(dat[[outlierY]], predict(quadratic))

    if(any(abs(rstandard(quadratic)) > STDRES) ){ #& abs(sd(residuals(quadratic))) > SDRES_MIN
      if(abs(sd(residuals(quadratic))) > SDRES_MIN){
        datOutQuad <- dat
        datOutQuad[[outlierY]][which(abs(rstandard(quadratic)) > STDRES)] <- NA
        datOutQuad[[outlierName]][which(abs(rstandard(quadratic)) > STDRES)] <- TRUE

        quadraticOut <- lm(get(outlierY) ~ poly(get(x), 2, raw = TRUE), data = datOutQuad)
        #RMSE
        quadraticRMSE <- Metrics::rmse(dat[[outlierY]], predict(quadratic))
        quadraticOutRMSE <- Metrics::rmse(na.exclude(datOutQuad[[outlierY]]), predict(quadraticOut))
        quadratic1 <- get(c("quadratic", "quadraticOut")[which(c(quadraticRMSE, quadraticOutRMSE) %in% min(quadraticRMSE, quadraticOutRMSE))])
        quadratic1.dat <-  datOutQuad

      } else{
        refslopemin <- coef(quadratic)[2]*100/3
        refslopemax <- coef(quadratic)[2]*100*3
        slopes <- sapply(1:(nrow(dat)-1), function(i) coef(lm(dat = dat[i:(i+1)], get(outlierY) ~ get(x)))[2]*100)
        if(slopes[1] > refslopemin){ slopes <- c(refslopemin*2, slopes)} else{slopes <- c(slopes[1],refslopemin*2, slopes[-1])}

        pos_slope <- c(which(abs(rstandard(quadratic)) > STDRES), which(abs(rstandard(quadratic)) > STDRES) + 1)
        pos_slope <- pos_slope[pos_slope <= nrow(dat)]


        if(any(slopes[pos_slope] < 0)){

          datOutQuad <- dat
          datOutQuad[[outlierY]][which(abs(rstandard(quadratic)) > STDRES)] <- NA
          datOutQuad[[outlierName]][which(abs(rstandard(quadratic)) > STDRES)] <- TRUE

          quadraticOut <- lm(get(outlierY) ~ poly(get(x), 2, raw = TRUE), data = datOutQuad)
          #RMSE
          quadraticRMSE <- Metrics::rmse(dat[[outlierY]], predict(quadratic))
          quadraticOutRMSE <- Metrics::rmse(na.exclude(datOutQuad[[outlierY]]), predict(quadraticOut))
          quadratic1 <- get(c("quadratic", "quadraticOut")[which(c(quadraticRMSE, quadraticOutRMSE) %in% min(quadraticRMSE, quadraticOutRMSE))])
          quadratic1.dat <-  datOutQuad


        } else{

          quadratic1 <- quadratic
          quadratic1.dat <- dat
        }

      }

    } else{
      quadratic1 <- quadratic
      quadratic1.dat <- dat
    }

    quadratic1RMSE <- quadraticRMSE #Metrics::rmse(quadratic1$model$`get(outlierY)`, predict(quadratic1))

    #cor.quadratic <- cor(dat[[tidyselect::all_of(y)]], predict(quadratic))
  } else{
    quadratic1 <- NA
    quadratic1RMSE <- NA
  }

  # select Model by using BIC
 # logistic1RMSE <- Metrics::rmse(logistic1$data$`get(tidyselect::all_of(y))`, predict(logistic1))
 # linear1RMSE <- Metrics::rmse(linear1$model$`get(tidyselect::all_of(y))`, predict(linear1))
 # quadratic1RMSE <- Metrics::rmse(quadratic1$model$`get(tidyselect::all_of(y))`, predict(quadratic1))
BIC_models <- BIC(logistic, linear, quadratic)

ModelName <- c("logistic1", "linear1", "quadratic1")[which(round(BIC_models$BIC,2) %in% min(round(BIC_models$BIC,2),na.rm = TRUE))]
#  ModelName <- c("logistic1", "linear1", "quadratic1")[which(round(c(logistic1RMSE, linear1RMSE, quadratic1RMSE),2) %in% min(round(c(logistic1RMSE, linear1RMSE, quadratic1RMSE),2),na.rm = TRUE))]
  if ("linear1" %in% ModelName) {ModelName = "linear1"} else if (all(c("logistic1", "quadratic1") %in% ModelName)) {ModelName = "logistic1"}  # if same correlation

  datName <- paste0(ModelName,".dat")

 #if(abs(sd(residuals(get(gsub(pattern = "1", x = ModelName, replacement = ""))))) > SDRES_MIN) dat[[outlierName]][which(abs(residuals(get(gsub(pattern = "1", x = ModelName, replacement = "")))/sd(residuals(get(gsub(pattern = "1", x = ModelName, replacement = ""))))) > STDRES)] <- TRUE
 #if(abs(sd(residuals(get(ModelName)))) > SDRES_MIN & any(abs(sd(residuals(get(ModelName)))/sd(residuals(get(ModelName)))) > STDRES)){
 #  dat[IDintern %in% get(get(ModelName)$call$data)[!is.na(get(y))][which(abs(residuals(get(ModelName))/sd(residuals(get(ModelName)))) > STDRES)]$IDintern][[outlierName]] <- TRUE
#}
  dat <- get(datName)
  dat <- data.table::setorder(dplyr::full_join(dat, dats[!IDintern %in% dat$IDintern], by = colnames(dats)), DilutionPoint)
  dat$color[dat[[outlierName]] %in% TRUE] <- "red"
  dat$Comment[dat[[outlierName]] %in% TRUE] <- paste0(dat$Comment[dat[[outlierName]] %in% TRUE], "_Outlier",abbr)
  dat$Comment[dat[[outlierName]] %in% FALSE] <- paste0(dat$Comment[dat[[outlierName]] %in% FALSE], "_NoOutlier",abbr)
  #dat[[outlierY]] <- dat[[y]]
  dat[[outlierY]][is.na(dat[, get(y)]) | dat[[outlierName]] %in% TRUE] <- NA



  modelNew <- if(ModelName %in% "logistic1"){drc::drm(get(outlierY) ~ get(x), fct = drc::L.3(), data = dat[color %in% "black",])
  } else if(ModelName %in% "linear1"){lm(get(outlierY) ~ get(x), data = dat[color %in% "black",])
      } else if(ModelName %in% "quadratic1"){lm(get(outlierY) ~ poly(get(x), 2, raw = TRUE), data = dat[color %in% "black",])}

  ## calculate correlation
  # cor.poly <- c(
  #   "cor.logistic" = cor.logistic,
  #   "cor.linear" = cor.linear,
  #   "cor.quadratic" = cor.quadratic
  # )
  # cor.max <-  names(which(cor.poly == max(cor.poly, na.rm = T)))
  #rm(list = c("logistic", "linear", "quadratic")[!c("logistic", "linear", "quadratic") %in% substr(cor.max,5, 100 )])



  Model <- ModelName

  Model = list("fit" = fitted(get(Model)), "coefficients" = coef(get(Model)), "std.residuals" = residuals(get(Model))/sd(residuals(get(Model))), "BIC" = BIC_models)



  tmp <- list(
    "tmp" = list(
      "model.name" = gsub(pattern = "1", x = ModelName, replacement = ""),
      "model" = Model,
      #"R2" = cor(dat[color %in% "black",][[y]], predict(modelNew))^2,
      #"aboveMinCor" = cor(dat[color %in% "black",][[y]], predict(modelNew))^2 > R2_MIN,
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

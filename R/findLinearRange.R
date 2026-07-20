#' Determine the linear dynamic range of a dilution or calibration series
#'
#' @description
#' `findLinearRange()` identifies the largest contiguous concentration or
#' dilution range that exhibits approximately linear behavior.
#'
#' The procedure first fits a linear regression model and iteratively removes
#' influential observations based on studentized residuals and Cook's distance.
#' Consecutive signals passing this filtering step are considered candidate
#' linear ranges.
#'
#' The final linear range must:
#'
#' * contain at least `min_feature` dilution points,
#' * include the midpoint of the response curve,
#' * show a positive monotonic response,
#' * form a continuous concentration/dilution interval.
#'
#' Once a suitable range has been identified, several performance metrics are
#' calculated, including:
#'
#' * slope (proportionality),
#' * adjusted R² (goodness-of-fit),
#' * Spearman correlation (monotonicity),
#' * maximum deviation from the fitted regression line (linearity).
#'
#' These metrics can subsequently be used to classify features as suitable or
#' unsuitable for quantitative analysis.
#'
#' @param dats Long-format data table containing a single dilution or
#'   concentration series.
#' @param x Character string specifying the transformed independent variable
#'   used for regression.
#' @param y Character string specifying the transformed dependent variable.
#' @param max_res Numeric threshold for studentized residuals. Observations with
#'   larger residuals are iteratively removed during range selection.
#' @param min_feature Minimum number of consecutive dilution points required for
#'   a valid linear range.
#' @param real_x Character string specifying the original (untransformed)
#'   concentration or dilution column.
#' @param slope_tol Numeric tolerance around an ideal slope of 1.
#'   Used for downstream linearity assessment.
#' @param delta_tol Maximum acceptable deviation (%) between observed and
#'   predicted values within the selected range.
#' @param rho_tol Allowed deviation from a perfect Spearman correlation
#'   coefficient of 1.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{SignalData}{Signal-level table containing the selected linear range and
#'   point-wise metrics.}
#'   \item{FeatureData}{Feature-level summary containing linearity metrics and range
#'   characteristics.}
#' }
#'
#' The signal-level output includes:
#'
#' \describe{
#'   \item{InRange}{Logical indicator showing whether a dilution point belongs
#'   to the selected linear range.}
#'   \item{Y_Range}{Response values restricted to the selected range.}
#'   \item{ResidualsInRange}{Residuals from the fitted linear model.}
#'   \item{delta}{Difference between observed and predicted values.}
#'   \item{delta_percent}{Relative deviation from the fitted model.}
#' }
#'
#' The feature-level output includes:
#'
#' \describe{
#'   \item{RangeStart}{First dilution point within the selected range.}
#'   \item{RangeEnd}{Last dilution point within the selected range.}
#'   \item{RangeLength}{Number of dilution points in the selected range.}
#'   \item{slope}{Slope of the fitted regression model.}
#'   \item{R2}{Adjusted R² of the fitted model.}
#'   \item{spearman_rho_inRange}{Monotonicity within the selected range.}
#'   \item{deltaMax_relative}{Maximum relative deviation from the fitted model.}
#' }
#'
#' @details
#' Unlike traditional calibration approaches that assess only the complete
#' concentration series, `findLinearRange()` identifies the portion of the
#' curve that behaves linearly. This allows robust evaluation of untargeted
#' metabolomics features whose dynamic range may be limited by detector
#' saturation, noise, or ion suppression.
#'
#' @seealso
#' \code{\link{create_output_findLinearRange}}
#'
#' @export


findLinearRange <- function(dats, x="DilutionPoint", y = "IntensityNorm",  max_res = 3, min_feature = 5, real_x, slope_tol = 0.2, delta_tol = 20, rho_tol = 0 ){




  dats$ResLR <- FALSE

  dat <- data.table::copy(dats)
  data.table::setorder(dat,DilutionPoint)
  dat <- dat[!is.na(get(y))]

  # calculate point nearest half max intensity
  int50 <- DescTools::Closest(x = dat[[y]] ,a = (min(dat[[y]]) + max(dat[[y]]))/2, which = TRUE, na.rm = T)
  if(length(int50) > 1) int50 <- max(int50)
  if(int50 == length(dat[[x]])) int50 <- length(dat[[x]]) -1
  if(int50 == 1) int50 = 2


  #create linear regression line going through int50
  #we <- rep(1, length(dat[[x]]))
  #we[(int50 - 1) : (int50 + 1)] <- 1000


  FIN <- FALSE

  while(FIN == FALSE) {

    we = rep(1, length(dat[[x]]))

    linearModel <- stats::lm(dat[[y]] ~ dat[[x]], weights = we)
    #quadratic <- lm(dat[[y]] ~ poly(dat[[x]], 2, raw = TRUE))

    fit <- stats::fitted(linearModel)

    ###use residuals

    n <- length(linearModel$fitted.values)
    p <- length(stats::coef(linearModel))


    threshold_cook <- 4 / (n - p - 1)


    fit_residuals <- round(stats::rstudent(linearModel), 1)
    cook <- round(stats::cooks.distance(linearModel), 2)

    lr <- !((abs(fit_residuals) > max_res & cook > threshold_cook))

    #dat$Residuals_weight_HalfmaxY = fit_residuals

    #lr <- abs(fit_residuals) < max_res

    if(all(lr %in% TRUE)) {
      FIN <- TRUE
      dat$ResLR <- TRUE

    } else if (all(lr %in% FALSE)) {
      FIN <- TRUE
    } else {
      dat <- dat[lr,]
    }

  }

  consNDX <- rle(lr)
  consNDX$position <- cumsum(consNDX$length)


  if(any(consNDX$length[which(consNDX$values %in% TRUE)] >= min_feature)){

    TRUEpos <- which( consNDX$position[consNDX$values %in% TRUE] >= int50  & consNDX$lengths[consNDX$values %in% TRUE] >= min_feature)

  } else{TRUEpos <- NULL}

  if(length(TRUEpos) > 0){


    if(length(TRUEpos) > 1) TRUEpos <- dplyr::last(TRUEpos)
    maxTrueRange <- (consNDX$position[consNDX$values %in% TRUE][TRUEpos] - consNDX$length[consNDX$values %in% TRUE][TRUEpos] +1) : consNDX$position[consNDX$values %in% TRUE][TRUEpos]
    maxTrueRange <- maxTrueRange[maxTrueRange!=0]

    dat$InRange <- dat$DilutionPoint >= dat$DilutionPoint[maxTrueRange[1]] & dat$DilutionPoint <= dat$DilutionPoint[utils::tail(maxTrueRange,1)]

    out <- create_output_findLinearRange(inRange = TRUE, data = dat, y = y, x = x, real_x = real_x, min_feature = min_feature)
    tmpGroup <- out[[1]]
    dat <- out[[2]]

    dat$Comment[dat$InRange %in% TRUE] <- unlist(apply(cbind(dat$Comment[dat$InRange %in% TRUE], "linearModel"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))

    # linear but not positive associated?
    if(any(dat$positiveSlope[dat$InRange %in% TRUE] %in% FALSE)){

      Range_TRUE <- which(dat$positiveSlope[dat$InRange %in% TRUE] %in% 1)
      Range_TRUE_list <- split(Range_TRUE, cumsum(c(1, diff(Range_TRUE) != 1)))
      Range_TRUE_list_Length <- lengths(Range_TRUE_list)

      indices <- 1:length( Range_TRUE_list)
      minsublist <- sapply(indices, function(i) min(Range_TRUE_list[[i]]))

      Range_TRUE_list <- lapply(indices, function(i) c(minsublist[i] -1,Range_TRUE_list[[i]]))
      Range_TRUE_list <- lapply(Range_TRUE_list, function(x) {x[x!=0]})

      if(length(Range_TRUE_list_Length[Range_TRUE_list_Length >= min_feature]) == 1 ){

        dat$InRange <- FALSE
        dat[unlist(Range_TRUE_list[which(Range_TRUE_list_Length >= min_feature)]), InRange := TRUE]
        #dat[unlist(Range_TRUE_list[which(Range_TRUE_list_Length < min_feature)]), InRange := FALSE]

        out <- create_output_findLinearRange(inRange = TRUE, data = dat, y = y, x = x, real_x = real_x, min_feature = min_feature)
        tmpGroup <- out[[1]]
        dat <- out[[2]]

      } else if(length(Range_TRUE_list_Length[Range_TRUE_list_Length >= min_feature]) > 1 &
                length(max(Range_TRUE_list_Length)) == 1){

        dat$InRange <- FALSE
        dat[unlist(Range_TRUE_list[which(Range_TRUE_list_Length == max(Range_TRUE_list_Length))]), InRange := TRUE]
        #dat[unlist(Range_TRUE_list[which(Range_TRUE_list_Length != max(Range_TRUE_list_Length))]), InRange := FALSE]

        out <- create_output_findLinearRange(inRange = TRUE, data = dat, y = y, x = x, real_x = real_x, min_feature = min_feature)
        tmpGroup <- out[[1]]
        dat <- out[[2]]

        tmpGroup$RangeFlag = "mutiple linear ranges"

      } else{


        tmpGroup$RangeFlag = "mutiple linear ranges, min 1 signal has negative slope"
      }}

  } else{

    dat$InRange <- FALSE
    out <- create_output_findLinearRange(inRange = FALSE, data = dat, y = y, x = x, real_x = real_x, min_feature = min_feature)
    tmpGroup <- out[[1]]
    dat <- out[[2]]

    dat$Comment <- unlist(apply(cbind(dat$Comment, "No linearModel"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))


  }

  dat <- data.table::setorder(dplyr::full_join(dat, dats[!IDintern %in% dat$IDintern], by = colnames(dats)), DilutionPoint)

  linearY <- "Y_Range"

  dat[[linearY]] <- dat[[y]]
  dat[[linearY]][is.na(dat[, get(y)]) | dat$InRange %in% FALSE | is.na(dat$InRange)] <- NA

  tmpGroup$Slope_within_Tolerance = abs(tmpGroup$slope - 1) <= slope_tol
  tmpGroup$Linearity_Criterion_Deviation = abs(tmpGroup$deltaMax_relative) <= delta_tol
  tmpGroup$Monotonicity_inRange =  abs(1 - round(tmpGroup$spearman_rho_inRange,2)) <= rho_tol
  tmpGroup$Monotonicity_complete = abs(1 - round(tmpGroup$spearman_rho_complete,2)) <= rho_tol


  dat$color = ifelse(dat$InRange %in% TRUE, "darkseagreen",dat$color)

  tmp <- list(dat, tmpGroup)
  return(tmp)
}



#' Calculate linearity metrics for a selected linear range
#'
#' @description
#' Internal helper function used by `findLinearRange()`.
#'
#' The function fits a linear regression model to the selected linear range and
#' calculates performance metrics describing monotonicity, proportionality,
#' goodness-of-fit, and deviation from linear behavior.
#'
#' @param inRange Logical value indicating whether a valid linear range was
#'   identified.
#' @param data Long-format data table containing a single dilution or
#'   concentration series.
#' @param y Character string specifying the transformed response variable.
#' @param x Character string specifying the transformed dilution or
#'   concentration variable.
#' @param real_x Character string specifying the original concentration or
#'   dilution column.
#' @param min_feature Minimum number of points required for a valid range.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{FeatureData}{Feature-level summary statistics.}
#'   \item{SignalData}{Signal-level table with fitted values, residuals, and
#'   deviation metrics.}
#' }
#'
#' Feature-level metrics include:
#'
#' \describe{
#'   \item{slope}{Slope of the fitted linear regression model.}
#'   \item{Intercept}{Regression intercept.}
#'   \item{R2}{Adjusted R².}
#'   \item{sigma}{Residual standard error.}
#'   \item{spearman_rho_inRange}{Spearman correlation within the selected
#'   range.}
#'   \item{spearman_rho_complete}{Spearman correlation across the complete
#'   series.}
#'   \item{deltaMax}{Maximum absolute deviation from the fitted regression
#'   line.}
#'   \item{deltaMax_relative}{Maximum relative deviation (%).}
#' }
#'
#' Signal-level metrics include:
#'
#' \describe{
#'   \item{predicted}{Predicted response values from the fitted model.}
#'   \item{ResidualsInRange}{Model residuals.}
#'   \item{delta}{Observed minus predicted values.}
#'   \item{delta_percent}{Relative deviation from predicted values.}
#'   \item{positiveSlope}{Indicator whether local response remains monotonic.}
#' }
#'
#' @details
#' This function is primarily intended for internal package use and is called
#' automatically by `findLinearRange()`.
#'
#' @keywords internal

create_output_findLinearRange <- function(inRange, data, y = y, x = x, real_x = real_x, min_feature = min_feature){



  if(inRange %in% TRUE){

    model <- stats::lm(get(y) ~ get(x), data = data[data$InRange %in% TRUE, ])
    fit <- stats::fitted(model)
    summary_model <- summary(model)

    spearman_rho <- stats::cor(data[[x]][data$InRange %in% TRUE], data[[y]][data$InRange %in% TRUE], method = "spearman")
    spearman_rho_all <- stats::cor(data[[x]], data[[y]], method = "spearman")

    y_obs <- exp(model$model[[1]])
    y_pred <- exp(stats::predict(model))

    rel_error <- abs((y_obs - y_pred) / y_obs)

    Deviation <- data[[y]][data$InRange %in% TRUE] - fit
    Deviation_perc <- abs(1 - exp(-Deviation))*100




    tmpGroup <- tibble::tibble(

      groupIndices = unique(data$groupIndices),
      RangeStart = data$DilutionPoint[data$InRange %in%TRUE][1],
      RangeStartY = data$Y[data$DilutionPoint %in% RangeStart],
      RangeStartX = data[[real_x]][data$DilutionPoint %in% RangeStart],
      RangeEnd = dplyr::last(data$DilutionPoint[data$InRange %in%TRUE]),
      RangeEndY = data$Y[data$DilutionPoint %in% RangeEnd],
      RangeEndX = data[[real_x]][data$DilutionPoint %in% RangeEnd],
      RangeLength = sum(data$InRange %in%TRUE),
      enoughPointsWithinRange = RangeLength >= min_feature,
      RangeFlag = NA,
      slope = stats::coef(model)[2],
      Intercept = stats::coef(model)[1],
      sigma = summary_model$sigma,
      R2 = summary_model$adj.r.squared,
      spearman_rho_inRange = spearman_rho,
      spearman_rho_complete = spearman_rho_all,
      positiveSlope_inRange = all(diff(data[[y]][data$InRange %in% TRUE]) > 0,na.rm = TRUE),
      positiveSlope_complete = all(diff(data[[y]]) > 0,na.rm = TRUE),
      deltaMax = max(abs(Deviation), na.rm = TRUE),
      deltaMax_relative = max(abs(Deviation_perc),na.rm = TRUE)


    )

    data <- data[, ':=' (
      positiveSlope = c(get(y)[1] < get(y)[2], diff(data[[y]]) > 0),
      R2 = ifelse(data$InRange %in% TRUE, tmpGroup$R2,NA),
      spearman_rho_linearRange = ifelse(data$InRange %in% TRUE, tmpGroup$spearman_rho_inRange,NA),
      spearman_rho = tmpGroup$spearman_rho_complete,
      deltaMax_relative = tmpGroup$deltaMax_relative




    )]

    data$predicted[data$InRange %in% TRUE] <- fit
    data$ResidualsInRange[data$InRange %in% TRUE] <- stats::resid(model)
    data$delta[data$InRange %in% TRUE] <- Deviation
    data$delta_percent[data$InRange %in% TRUE] <- Deviation_perc


  } else {

    spearman_rho_all <- stats::cor(data[[x]], data[[y]], method = "spearman")

    tmpGroup <- tibble::tibble(

      groupIndices = unique(data$groupIndices),
      RangeStart = NA,
      RangeStartY = NA,
      RangeStartX = NA,
      RangeEnd = NA,
      RangeEndY = NA,
      RangeEndX = NA,
      RangeLength = NA,
      enoughPointsWithinRange = NA,
      RangeFlag = NA,
      slope = NA,
      Intercept = NA,
      sigma = NA,
      R2 = NA,
      spearman_rho_inRange = NA,
      spearman_rho_complete = spearman_rho_all,
      positiveSlope_inRange = NA,
      positiveSlope_complete = all(diff(data[[y]]) > 0,na.rm = TRUE),
      #monoton = NA,
      deltaMax = NA,
      deltaMax_relative = NA


    )

    data <- data[, ':=' (
      positiveSlope = NA,
      R2 = NA,
      spearman_rho_linearRange = NA,
      spearman_rho = tmpGroup$spearman_rho_complete,
      predicted = NA,
      ResidualsInRange = NA,
      delta = NA,
      delta_percent = NA,
      deltaMax_relative = NA

    )]



  }

  return(list(tmpGroup, data))

}
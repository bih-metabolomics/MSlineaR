#' determine metrics for linearity and linear range of dilution/ calibration curve
#'
#' @description `findLinearRange()` first determines a suitable range by forcing
#' a linear regression model through halfmax intensity, calculating the residuals and
#' removing all dilution points where the residuals are bigger than predefined value `max_res`.
#' Using this range the method calculates metrics for monotonicity using spearman rho,
#' linearity by using deviation of observed to predicted values from linear regression model,
#' proportionality across range by calculating the slope of linear regression model
#' and the goodness-of-fit by calculating the adjusted R^2.
#' @param x String; Column name of the log transformed independent variable
#' @param dats Long format data frame or data table for one metabolite including information about dependent and independent variable.
#' Further necessary columns are: groupIndices, DilutionPoint (1:x), IDintern, color.
#' @param max_res Integer; points of serial diluted/concentrated series,
#' which are less than `max_res` are considered for checking of linearity.
#' Default to 3.
#' @param min_feature Integer, ranging between 3 and maximum number of dilutions/concentrations.
#' Minimum number of points present in one serial diluted/concentrated series
#' marked as linear to consider this signal as linear. Default to 6, according to EMA guidelines2022.
#' @param real_x String; Column name of the untransformed independent variable
#' @param slope_tol Integer, allowed tolerance for slope deviation. Default to 0.15
#' @param delta_tol Integer, allowed tolerance for maximum deviation of observed
#' values to predicted value from linear regression model. Default set to 0.182 (~ 20%).
#' @param y String; Column name of the log transformed dependent variable
#'
#' @return list of two data frames. One with results per dilution point, one with result of features.
#'
#'
#' @export
#'

#' @import dplyr
#' @import data.table
#' @importFrom tibble tibble
#' @importFrom Matrix tail
#' @importFrom DescTools Closest
#' @importFrom stats fitted lm residuals

findLinearRange <- function(dats, x="DilutionPoint", y = "IntensityNorm",  max_res = 3, min_feature = 5, real_x, slope_tol = 0.2, delta_tol = 20, rho_tol = 0 ){

  create_output_findLinearRange <- function(inRange, data, y = y, x = x, real_x = real_x, min_feature = min_feature){



    if(inRange %in% TRUE){

      model <- lm(get(y) ~ get(x), data = data[data$InRange %in% TRUE, ])
      fit <- fitted(model)
      summary_model <- summary(model)

      spearman_rho <- cor(data[[x]][data$InRange %in% TRUE], data[[y]][data$InRange %in% TRUE], method = "spearman")
      spearman_rho_all <- cor(data[[x]], data[[y]], method = "spearman")

      Deviation <- data[[y]][data$InRange %in% TRUE] - fit
      Deviation_perc <- (exp(Deviation) -1)*100




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
        slope = coef(model)[2],
        Intercept = coef(model)[1],
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
      data$ResidualsInRange[data$InRange %in% TRUE] <- resid(model)
      data$delta[data$InRange %in% TRUE] <- Deviation
      data$delta_percent[data$InRange %in% TRUE] <- Deviation_perc


    } else{

      spearman_rho_all <- cor(data[[x]], data[[y]], method = "spearman")

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


  FIN = FALSE

  while(FIN == FALSE) {

    we = rep(1, length(dat[[x]]))

    linearModel <- lm(dat[[y]] ~ dat[[x]] , weights = we)
    #quadratic <- lm(dat[[y]] ~ poly(dat[[x]], 2, raw = TRUE))

    fit <- fitted(linearModel)

    ###use residuals

    n <- length(linearModel$fitted.values)
    p <- length(coef(linearModel))


    threshold_cook <- 4 / (n - p - 1)


    fit_residuals <- round(rstudent(linearModel),1)
    cook <- round(cooks.distance(linearModel),2)

    lr <- !((abs(fit_residuals) > max_res & cook > threshold_cook))

    #dat$Residuals_weight_HalfmaxY = fit_residuals

    #lr <- abs(fit_residuals) < max_res

    if(all(lr %in% TRUE)) {
      FIN = TRUE
      dat$ResLR <- TRUE

    } else if (all(lr %in% FALSE)) {
      FIN = TRUE
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

    dat$InRange <- dat$DilutionPoint >= dat$DilutionPoint[maxTrueRange[1]] & dat$DilutionPoint <= dat$DilutionPoint[tail(maxTrueRange,1)]

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



#' calculate linear regression from choosen range and generate output table for function `findLinearRange()`
#'
#' @param inRange Boolean; Does the feature has a choosen range?
#' @param data Long format data frame or data table for one metabolite generated by function `findLinearRange()`
#' @param y String; Column name of the log transformed dependent variable
#' @param x String; Column name of the log transformed independent variable
#'
#' @return list including two tables with metrics for each signal in the first one and information about the whole feature in the second one

#'

# create_output_findLinearRange <- function(inRange, data, y = y, x = x, real_x = real_x, min_feature = min_feature){
#
#
#
#   if(inRange %in% TRUE){
#
#     model <- lm(get(y) ~ get(x), data = data[data$InRange %in% TRUE, ])
#     fit <- fitted(model)
#     summary_model <- summary(model)
#
#     spearman_rho <- cor(data[[x]][data$InRange %in% TRUE], data[[y]][data$InRange %in% TRUE], method = "spearman")
#
#     Deviation <- data[[y]][data$InRange %in% TRUE] - fit
#     Deviation_perc <- (exp(Deviation) -1)*100
#
#
#
#
#     tmpGroup <- tibble::tibble(
#
#       groupIndices = unique(data$groupIndices),
#       RangeStart = data$DilutionPoint[data$InRange %in%TRUE][1],
#       RangeStartY = data$Y[data$DilutionPoint %in% RangeStart],
#       RangeStartX = data[[real_x]][data$DilutionPoint %in% RangeStart],
#       RangeEnd = last(data$DilutionPoint[data$InRange %in%TRUE]),
#       RangeEndY = data$Y[data$DilutionPoint %in% RangeEnd],
#       RangeEndX = data[[real_x]][data$DilutionPoint %in% RangeEnd],
#       RangeLength = sum(data$InRange %in%TRUE),
#       enoughPointsWithinRange = RangeLength >= min_feature,
#       RangeFlag = NA,
#       slope = coef(model)[2],
#       Intercept = coef(model)[1],
#       R2 = summary_model$adj.r.squared,
#       spearman_rho = spearman_rho,
#       positiveSlope = all(diff(data[[y]][data$InRange %in% TRUE]) > 0,na.rm = TRUE),
#       deltaMax = max(abs(Deviation), na.rm = TRUE),
#       deltaMax_relative = max(abs(Deviation_perc),na.rm = TRUE)
#
#
#     )
#
#     data <- data[, ':=' (
#       positiveSlope = c(get(y)[1] < get(y)[2], diff(data[[y]]) > 0),
#       R2 = ifelse(data$InRange %in% TRUE, tmpGroup$R2,NA)
#
#
#
#
#     )]
#
#     data$predicted[data$InRange %in% TRUE] <- fit
#     data$ResidualsInRange[data$InRange %in% TRUE] <- resid(model)
#     data$delta[data$InRange %in% TRUE] <- Deviation
#     data$delta_percent[data$InRange %in% TRUE] <- Deviation_perc
#
#
#   } else{
#
#     tmpGroup <- tibble::tibble(
#
#       groupIndices = unique(data$groupIndices),
#       RangeStart = NA,
#       RangeStartY = NA,
#       RangeStartX = NA,
#       RangeEnd = NA,
#       RangeEndY = NA,
#       RangeEndX = NA,
#       RangeLength = NA,
#       enoughPointsWithinRange = NA,
#       RangeFlag = NA,
#       slope = NA,
#       Intercept = NA,
#       R2 = NA,
#       spearman_rho = NA,
#       monoton = NA,
#       deltaMax = NA,
#       deltaMax_relative = NA
#
#
#     )
#
#     data <- data[, ':=' (
#       positiveSlope = NA,
#       R2 = NA,
#       predicted = NA,
#       ResidualsInRange = NA,
#       deta = NA,
#       delta_percent = NA
#
#     )]
#
#
#
#   }
#
#   return(list(tmpGroup, data))
#
# }
















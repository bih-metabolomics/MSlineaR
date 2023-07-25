#' Title
#'
#' @param dat
#' @param x
#' @param y
#' @param res
#'
#' @return
#' @export
#'
#' @examples
findLinearRange <- function(dats, x="DilutionPoint", y = "IntensityNorm",  sd_res_factor = 2, min_feature, real_x){#modelObject
  #browser()

  dat <- data.table::copy(dats)
  data.table::setorder(dat,DilutionPoint)
  dat <- dat[!is.na(get(y))]
  #modelObject <- unique(dat[[modelObject]])
  #modelObject <- unlist(modelObject, recursive = F)
  #int50 <- DescTools::Closest(x = dat[[y]] ,a = (max(modelObject$fit) -min(modelObject$fit))/2 + min(modelObject$fit), which = TRUE, na.rm = T)
  int50 <- DescTools::Closest(x = dat[[y]] ,a = (max(dat[[y]]) -min(dat[[y]]))/2 + min(dat[[y]]), which = TRUE, na.rm = T)

  if(length(int50) > 1) int50 <- max(int50)
  if(int50 == length(dat[[x]])) int50 <- length(dat[[x]]) -1
  if(int50 == 1) int50 = 2
  #dat$color[int50] <-  "green"



  #create linear regression line going through int50
  we <- rep(1, length(dat[[x]]))
  we[(int50 - 1) : (int50 + 1)] <- 10
  #we = NULL


  linearRange <- lm(dat[[y]] ~ dat[[x]], weights = we)
  ablineIntensity <- fitted(linearRange)

  ###use residuals

  std_residuals <- rstandard(linearRange)
  sd_residuals <- abs(sd_res_factor*sd(std_residuals[which(abs(std_residuals) < 3)]))
  if(sd_residuals < 1) sd_residuals <- 1


  ##use cooks distance
  cook <- cooks.distance(linearRange)
  cookref <- dat$DilutionPoint[which(cook > 1)]

  lr <- !(abs(std_residuals) >= ceiling(sd_residuals*10)/10 & dat$DilutionPoint %in% cookref)

  #lr <- !(abs(std_residuals) >= ceiling(sd_residuals*10)/10 & dat$DilutionPoint %in% cookref)

   consNDX <- rle(lr)

  consNDX$position <- cumsum(consNDX$length)


  if(any(consNDX$length[which(consNDX$values %in% TRUE)]>= min_feature)){ #& any(consNDX$position[consNDX$values %in% TRUE] >= int50)

    TRUEpos <- Position(function(fi) fi >= int50, consNDX$position[consNDX$values %in% TRUE], right = TRUE)
    maxTrueRange <- (consNDX$position[consNDX$values %in% TRUE
    ][TRUEpos] - consNDX$length[consNDX$values %in% TRUE][TRUEpos] +1) : consNDX$position[consNDX$values %in% TRUE][TRUEpos]
    maxTrueRange <- maxTrueRange[maxTrueRange!=0]
    # if(length(consNDX$length[which(consNDX$values %in% TRUE)]))


    tmpGroup <- tibble::tibble(
      groupIndices = unique(dat$groupIndices),
      linear = TRUE,
      LRStart = dat$DilutionPoint[maxTrueRange[1]],
      LRStartY = dat$Y[maxTrueRange[1]],
      LRStartX = dat[[real_x]][maxTrueRange[1]],
      LREnd = dat$DilutionPoint[tail(maxTrueRange,1)],
      LREndY = dat$Y[tail(maxTrueRange,1)],
      LREndX = dat[[real_x]][tail(maxTrueRange,1)],
      LRLength = length(maxTrueRange),
      enoughPointsWithinLR = LRLength >= min_feature,
      LRFlag = NA

    )

    dat[, ':=' (IsLinear = DilutionPoint >= tmpGroup$LRStart & DilutionPoint <= tmpGroup$LREnd,
                IsPositivAssociated = c(get(y)[1] < get(y)[2], (get(y)[-1] - data.table::shift(get(y), 1, type = "lag")[-1]) > 0),
                #modelFit = modelObject$fit,
                #ablineLimit1 = limitup,
                #ablineLimit2 = limitdown,
                Residuals = std_residuals

    )]
    dat$color[dat$IsLinear %in% TRUE] <- "darkseagreen"
    dat$R2[dat$IsLinear %in% TRUE] <- summary(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))$r.squared
    dat$abline[dat$IsLinear %in% TRUE] = fitted(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))
    tmpGroup$Intercept <- lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ])$coefficients[[1]]
    tmpGroup$slope <- lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ])$coefficients[[2]]
    tmpGroup$R2 <- summary(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))$r.squared

    dat$Comment[dat$IsLinear %in% TRUE] <- unlist(apply(cbind(dat$Comment[dat$IsLinear %in% TRUE], "linearRange"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))

    # linear but not positive associated?
    if(any(dat$IsPositivAssociated[tmpGroup$LRStart : tmpGroup$LREnd] %in% FALSE)){

      LR_TRUE <- which(dat$IsPositivAssociated %in% 1)
      LR_TRUE_list <- split(LR_TRUE, cumsum(c(1, diff(LR_TRUE) != 1)))
      LR_TRUE_list_Length <- lengths(LR_TRUE_list)

      indices <- 1:length( LR_TRUE_list)
      minsublist <- sapply(indices, function(i) min(LR_TRUE_list[[i]]))

      LR_TRUE_list <- lapply(indices, function(i) c(minsublist[i] -1,LR_TRUE_list[[i]]))
      LR_TRUE_list <- lapply(LR_TRUE_list, function(x) {x[x!=0]})

      if(any(lengths(LR_TRUE_list) >= min_feature & length(LR_TRUE_list) == 1)){


        dat[unlist(LR_TRUE_list[which(LR_TRUE_list_Length >= min_feature)]), IsLinear := TRUE]
        dat[unlist(LR_TRUE_list[which(LR_TRUE_list_Length < min_feature)]), IsLinear := FALSE]
        dat$color[dat$IsLinear %in% TRUE] <- "darkseagreen"
        dat$color[dat$IsLinear %in% FALSE] <- "black"
        dat$R2[dat$IsLinear %in% TRUE] <- summary(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))$r.squared
        dat$abline[dat$IsLinear %in% TRUE] = fitted(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))

        tmpGroup$LRStart = dat$DilutionPoint[dat$IsLinear %in% TRUE][1]
        tmpGroup$LRStartY = dat$Y[dat$IsLinear %in% TRUE][1]
        tmpGroup$LRStartX =  dat[[real_x]][dat$IsLinear %in% TRUE][1]
        tmpGroup$LREnd = data.table::last(dat$DilutionPoint[dat$IsLinear %in% TRUE])
        tmpGroup$LREndY = data.table::last(dat$Y[dat$IsLinear %in% TRUE])
        tmpGroup$LREndX = data.table::last(dat[[real_x]][dat$IsLinear %in% TRUE])
        tmpGroup$LRLength = sum(dat$IsLinear)
        tmpGroup$enoughPointsWithinLR = tmpGroup$LRLength >= min_feature
        tmpGroup$Intercept <- lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ])$coefficients[[1]]
        tmpGroup$slope <- lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ])$coefficients[[2]]
        tmpGroup$R2 <- summary(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))$r.squared

        } else if(any(lengths(LR_TRUE_list) >= min_feature & length(LR_TRUE_list) > 1 & length(max(LR_TRUE_list_Length)) == 1)){

          dat[unlist(LR_TRUE_list[which(LR_TRUE_list_Length == max(LR_TRUE_list_Length))]), IsLinear := TRUE]
          dat[unlist(LR_TRUE_list[which(LR_TRUE_list_Length != max(LR_TRUE_list_Length))]), IsLinear := FALSE]
          dat$color[dat$IsLinear %in% TRUE] <- "darkseagreen"
          dat$color[dat$IsLinear %in% FALSE] <- "black"
          dat$R2[dat$IsLinear %in% TRUE] <- summary(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))$r.squared
          dat$abline[dat$IsLinear %in% TRUE] = fitted(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))


          tmpGroup$LRStart = dat$DilutionPoint[dat$IsLinear %in% TRUE][1]
          tmpGroup$LRStartY = dat$Y[dat$IsLinear %in% TRUE][1]
          tmpGroup$LRStartX =  dat[[real_x]][dat$IsLinear %in% TRUE][1]
          tmpGroup$LREnd = data.table::last(dat$DilutionPoint[dat$IsLinear %in% TRUE])
          tmpGroup$LREndY = data.table::last(dat$Y[dat$IsLinear %in% TRUE])
          tmpGroup$LREndX = data.table::last(dat[[real_x]][dat$IsLinear %in% TRUE])
          tmpGroup$LRLength = sum(dat$IsLinear)
          tmpGroup$enoughPointsWithinLR = tmpGroup$LRLength >= min_feature
          tmpGroup$Intercept <- lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ])$coefficients[[1]]
          tmpGroup$slope <- lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ])$coefficients[[2]]
          tmpGroup$R2 <- summary(lm(get(y) ~ get(x), data = dat[color %in% "darkseagreen", ]))$r.squared
          tmpGroup$LRFlag = "mutiple linear ranges"

        } else{

          dat$IsLinear = FALSE
          dat$color[dat$IsLinear %in% TRUE] <- "darkseagreen"
          dat$color[dat$IsLinear %in% FALSE] <- "black"
          dat$R2 <- NA
          dat$abline <- NA

          tmpGroup$linear = FALSE
          tmpGroup$LRStart = NA
          tmpGroup$LRStartY = NA
          tmpGroup$LRStartX =  NA
          tmpGroup$LREnd = NA
          tmpGroup$LREndY = NA
          tmpGroup$LREndX = NA
          tmpGroup$LRLength = NA
          tmpGroup$enoughPointsWithinLR = NA
          tmpGroup$LRFlag = NA
          tmpGroup$Intercept <- NA
          tmpGroup$slope <- NA
          tmpGroup$R2 <- NA
      }}

      } else{

        tmpGroup <- tibble::tibble(
          groupIndices = unique(dat$groupIndices),
          linear = FALSE,
          LRStart = NA,
          LRStartY = NA,
          LRStartX = NA,
          LREnd = NA,
          LREndY = NA,
          LREndX = NA,
          LRLength = NA,
          enoughPointsWithinLR = NA,
          LRFlag = NA,
          Intercept = NA,
          slope = NA,
          R2 = NA

        )

        dat[, ':=' (IsLinear = FALSE,
                    IsPositivAssociated = c(get(y)[1] < get(y)[2], (get(y)[-1] - data.table::shift(get(y), 1, type = "lag")[-1]) > 0),
                    #modelFit = modelObject$fit,
                    #ablineLimit1 = limitup,
                    #ablineLimit2 = limitdown,
                    Residuals = std_residuals,
                    R2 = NA,
                    abline = NA
        )]
        dat$Comment[dat$IsLinear %in% FALSE] <- unlist(apply(cbind(dat$Comment[dat$IsLinear %in% FALSE], "NolinearRange"), 1, function(x) paste(x[!is.na(x)], collapse = "_")))

      }
  dat <- data.table::setorder(dplyr::full_join(dat, dats[!IDintern %in% dat$IDintern], by = colnames(dats)), DilutionPoint)

  linearY <- "Y_LR"
  dat[[linearY]] <- dat[[y]]
  dat[[linearY]][is.na(dat[, get(y)]) | dat$IsLinear %in% FALSE | is.na(dat$IsLinear)] <- NA
  dat$Intercept <- tmpGroup$Intercept
  dat$slope <- tmpGroup$slope
  dat$LRStartY <- tmpGroup$LRStartY
  dat$LRStart <- tmpGroup$LRStart
  dat$LREndY <- tmpGroup$LREndY
  dat$LREnd <- tmpGroup$LREnd



#dat <- subset(dat, select = -fittingModel)
  tmp <- list(dat, tmpGroup)
  return(tmp)
}










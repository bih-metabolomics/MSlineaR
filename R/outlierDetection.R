#' Title
#'
#' @param ID
#' @param modelObject
#' @param res
#' @param count
#'
#' @return
#' @export
#'
#' @examples
outlier <- function(ID, modelObject, res, count){

  residual.all <- abs(residuals(modelObject)/sd(residuals(modelObject)))
  residual <- sort(residual.all,decreasing = T)
  if ( count > length(residual.all)) count = length(residual.all)
  residual <- residual[1:count]
  residual <- residual[residual > res]

  if(length(residual) > 0) {
    rs <- ldply(residual, function(x){
      data.frame(
        "ID" = ID,
        "outlier" = TRUE,
        "outlier_position" = which(residual.all == x)
      )
    })
    #   data.frame(
    #   "ID" = ID,
    #   "outlier" = TRUE,
    #   "outlier_position" = which(residual == max(residual))
    # )
  } else {
    rs <- data.frame(
      "ID" = ID,
      "outlier" = FALSE,
      "outlier_position" = NA
    )
  }

  return(rs)
}

#' Title
#'
#' @param dat
#' @param y
#' @param x
#' @param model
#' @param res
#' @param threshCor
#' @param numboutlier
#'
#' @return
#' @export
#'
#' @examples
outlierDetection <- function(dat, y="IntensityNorm", x="DP", model=c("logistic", "linear", "quadratic"), res=2, threshCor=0.99, numboutlier = 1){
  #browser()
  dat <- dat %>% arrange(DilutionPoint)
  dataOutlier <- dat %>% drop_na(y)
  dataOutlier$Outlier <- NA
  bestModel <- chooseModel(dat, y, x, model)

  dataModel <- data.frame(
    "ID" = names(bestModel),
    "ModelFit" = map(bestModel, 1) %>%  unlist(use.names = F),
    "correlation" = map(bestModel, 3) %>%  unlist(use.names = F)
  )

  if (dataModel$correlation < threshCor) {
    outlierSum <- outlier(ID = unique(dat$ID),modelObject =  bestModel[[1]][[2]], res, count = numboutlier)
    if (any(outlierSum$outlier %in% TRUE)) {
      dataOutlier$color[outlierSum$outlier_position] <- "red"
      dataOutlier$pch[outlierSum$outlier_position] <- 19
      dataOutlier$Outlier[outlierSum$outlier_position] <- TRUE
    }
  } else{
    outlierSum <- outlier(ID = unique(dat$ID),modelObject =  bestModel[[1]][[2]], res = 100, count = numboutlier)
  }

  dataOutlier <- full_join(dataOutlier,
                           dat %>% dplyr::select(-color,-pch),
                           by = names(dat)[-which(arr.ind = T,names(dat) %in% c("color", "pch" ))])
  tmp <- list("tmp" = list("Model" = dataModel, "Outlier" = outlierSum, "dataOutlier" = dataOutlier))
  names(tmp) <- unique(dat$ID)
  return(tmp)

}

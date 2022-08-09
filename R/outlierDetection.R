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
outlier <- function(dat, modelObject, res, count){

  residual.all <- abs(residuals(modelObject)/sd(residuals(modelObject)))
  residual <- sort(residual.all,decreasing = T)
  if ( count > length(residual.all)) count = length(residual.all)
  residual <- residual[1:count]
  residual <- residual[residual > res]


  dat$outlier = residual.all %in% residual
  #"outlier_position" = which(residual.all %in% residual),
  dat$residuals = residual.all

  return(dat)
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
outlierDetection <- function(dat, y="IntensityNorm", x="DilutionPoint", model=c("logistic", "linear", "quadratic"), res=2, threshCor=0.99, numboutlier = 1){
  #browser()
  dat <- dat %>% arrange(DilutionPoint)
  dataOutlier <- dat %>% drop_na(all_of(y))
  #dataOutlier$Outlier <- NA
  bestModel <- chooseModel(dat, all_of(y), x, model)

  dataModel <- tibble(
    "groupIndices" = as.integer(names(bestModel)),
    #"ID" = names(bestModel),
    "ModelFit" = map(bestModel, 1) %>%  unlist(use.names = F),
    "correlationModel" = map(bestModel, 3) %>%  unlist(use.names = F)
  )

  if (dataModel$correlationModel < threshCor) {
    dataOutlier <- outlier(dataOutlier,modelObject =  bestModel[[1]][[2]], res, count = numboutlier)
    if (any(dataOutlier$outlier %in% TRUE)) {
      dataOutlier$color[dataOutlier$outlier %in% TRUE] <- "red"
      dataOutlier$pch[dataOutlier$outlier %in% TRUE] <- 19
      dataOutlier$Comment[dataOutlier$outlier %in% TRUE] <- "outlier"
    }
  } else{
    dataOutlier <- outlier(dataOutlier , modelObject =  bestModel[[1]][[2]], res = 100, count = numboutlier)
  }

  Outlier <- full_join(dataOutlier, dataModel, by = "groupIndices") |> dplyr::select(IDintern, groupIndices,Comment,outlier,pch, color, residuals, ModelFit, correlationModel )


  return(Outlier)

}

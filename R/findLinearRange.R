#' Title
#'
#' @param dat
#' @param x
#' @param y
#' @param modelObject
#' @param res
#'
#' @return
#' @export
#'
#' @examples
findLinearRange <- function(dat, x="DilutionPoint", y = "IntensityNorm", modelObject, res = 0.2, minConsecutives){
  #browser()

  modelObject <- unique(dat[[modelObject]])
  int50 <- DescTools::Closest(x = dat[[y]] ,a = max( fitted(modelObject[[1]]))/2, which = TRUE, na.rm = T)

  dat$color[int50] <-  "green"
  dat$pch[int50] <- 19

  #create linear regression line going through int50
  linearRange <- lm(fitted(modelObject[[1]])[(int50 - 1) : (int50 + 1)] ~ dat[[x]][(int50 - 1) : (int50 + 1)])
  ablineIntensity <- coef(linearRange)[1] + coef(linearRange)[2]*dat[[x]]

  ndx <- which(abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))) < res)

  #if(length(ndx) >= 2){
    tmp <- tibble(
      IDintern = dat$IDintern,
      groupIndices = dat$groupIndices,
      linearRangeStart = dat$DilutionPoint[ndx[1]],
      linearRangeEnd = dat$DilutionPoint[tail(ndx,1)],
      IslinearRange = dat$DilutionPoint >= linearRangeStart & dat$DilutionPoint <= linearRangeEnd,
      IsPositivAssociated = (dat$IntensityNorm - lag(dat$IntensityNorm)) > 0,
      linearRange = linearRangeEnd - linearRangeStart + 1,
      enoughPointsWithinLinearRange = linearRange >= minConsecutives,
      color = dat$color,
      pch = dat$pch,
      modelFit = fitted(modelObject[[1]]),
      ablineFit = ablineIntensity,
      ablineRes = abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))),
      Comment = "linearRange found"

    )
#
#   } else{
#     tmp <- tibble(
#       IDintern = dat$IDintern,
#       groupIndices = dat$groupIndices,
#       linearRangeStart = NA,
#       linearRangeEnd = NA,
#       IslinearRange = FALSE,
#       IsPositivAssociated = (dat$IntensityNorm - lag(dat$IntensityNorm)) > 0,
#       linearRange = NA,
#       enoughPointsWithinLinearRange = NA,
#       color = dat$color,
#       pch = dat$pch,
#       modelFit = fitted(modelObject[[1]]),
#       ablineFit = ablineIntensity,
#       ablineRes = abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))),
#       Comment = "no linearRange found"
#
#     )
#     }

  return(tmp)
}

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

  setorder(dat,DilutionPoint)
  dat <- dat[color %in% "black"]
  modelObject <- unique(dat[[modelObject]])
  int50 <- DescTools::Closest(x = dat[[y]] ,a = max( fitted(modelObject[[1]]))/2, which = TRUE, na.rm = T)

  dat$color[int50] <-  "green"
  dat$pch[int50] <- 19

  #create linear regression line going through int50
  linearRange <- lm(fitted(modelObject[[1]])[(int50 - 1) : (int50 + 1)] ~ dat[[x]][(int50 - 1) : (int50 + 1)])
  ablineIntensity <- coef(linearRange)[1] + coef(linearRange)[2]*dat[[x]]

  #ndx <- which(abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))) < res)
  consNDX <- rle(abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))) < res)

  consNDX$position <- cumsum(consNDX$length)

  if(any(consNDX$length[which(consNDX$values %in% TRUE)]>= 3)){

  TRUEpos <- Position(function(fi) fi >= int50, consNDX$position[consNDX$values %in% TRUE], right = TRUE)
  maxTrueRange <- (consNDX$position[consNDX$values %in% TRUE
                                   ][TRUEpos] - consNDX$length[consNDX$values %in% TRUE][TRUEpos]) : consNDX$position[consNDX$values %in% TRUE][TRUEpos]
maxTrueRange <- maxTrueRange[maxTrueRange!=0]
 # if(length(consNDX$length[which(consNDX$values %in% TRUE)]))



    tmp <- tibble(
      IDintern = dat$IDintern,
      groupIndices = dat$groupIndices,
      linearRangeStart = dat$DilutionPoint[maxTrueRange[1]],
      linearRangeEnd = dat$DilutionPoint[tail(maxTrueRange,1)],
      IslinearRange = dat$DilutionPoint >= linearRangeStart & dat$DilutionPoint <= linearRangeEnd,
      IsPositivAssociated = (dat$IntensityNorm - lag(dat$IntensityNorm)) > 0,
      linearRange = length(maxTrueRange),
      enoughPointsWithinLinearRange = linearRange >= minConsecutives,
      color = dat$color,
      pch = dat$pch,
      modelFit = fitted(modelObject[[1]]),
      ablineFit = ablineIntensity,
      ablineRes = abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))),
      Comment = "linearRange found"

    )

  } else{
    tmp <- tibble(
      IDintern = dat$IDintern,
      groupIndices = dat$groupIndices,
      linearRangeStart = NA,
      linearRangeEnd = NA,
      IslinearRange = FALSE,
      IsPositivAssociated = (dat$IntensityNorm - lag(dat$IntensityNorm)) > 0,
      linearRange = NA,
      enoughPointsWithinLinearRange = NA,
      color = dat$color,
      pch = dat$pch,
      modelFit = fitted(modelObject[[1]]),
      ablineFit = ablineIntensity,
      ablineRes = abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))),
      Comment = "no linearRange found"

    )
    }

  return(tmp)
}

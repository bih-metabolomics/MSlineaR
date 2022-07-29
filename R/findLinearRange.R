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
findLinearRange <- function(dat, x="DP", y = "IntensityNorm", modelObject, res = 0.1){
  #browser()

  int50 <- DescTools::Closest(x = dat[[y]] ,a = max( fitted(modelObject))/2, which = TRUE)

  dat$color[int50] <-  "green"
  dat$pch[int50] <- 19

  #create linear regression line going through int50
  linearRange <- lm(fitted(modelObject)[(int50 - 1) : (int50 + 1)] ~ dat[[x]][(int50 - 1) : (int50 + 1)])
  ablineIntensity <- coef(linearRange)[1] + coef(linearRange)[2]*dat[[x]]

  ndx <- which(abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))) < res)

  if(length(ndx) >= 2){
    dat$linearRangeStart <- dat$DilutionPoint[ndx[1]]
    dat$linearRangeEnd <- dat$DilutionPoint[tail(ndx,1)]
    dat$linearRange[dat$DilutionPoint >= dat$linearRangeStart & dat$DilutionPoint <= dat$linearRangeEnd] <- TRUE
    dat$color[int50] <- "green"
    dat$fitted <- fitted(modelObject)
    dat$fittedLM <- ablineIntensity
  } else{
    dat$linearRangeStart <- NA
    dat$linearRangeEnd <- NA
    dat$linearRange <- NA
    dat$fitted <- NA
    dat$fittedLM <- NA
    dat$Comment = "no linear Range found"}

  return(dat)
}

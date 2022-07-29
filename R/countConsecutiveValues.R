#' Title
#'
#' @param dats
#' @param minConsecutives
#' @param y
#'
#' @return
#' @export
#'
#' @examples
consecutiveVali <- function(dats, minConsecutives=5, y="IntensityNorm"){

  dat <- dats %>% arrange(DilutionPoint)
  if(any(is.na(dat[[y]]))){

    #remove NAs at start and end
    naLimits <- rle(is.na(dat[[y]]))
    if(naLimits$values[1] %in% TRUE){
      dat <- dat[-c(1:naLimits$lengths[1]), ]
    }

    if(naLimits$values[length(naLimits)] %in% TRUE){
      dat <- dat[-c((nrow(dat) - naLimits$lengths[length(naLimits)] + 1) : nrow(dat)), ]
    }
  }


  valid.dat <- dat %>% filter(! Comment %in% c("bigger", "smaller", "last", "first"))
  cons <- rle(diff(valid.dat$DilutionPoint) == 1)

  tmp <- data.frame("ID" = unique(dat$ID), "enoughPoints" = any(cons$lengths >= (minConsecutives-1) & cons$values == TRUE & sum(dat$color %in% "black") >= minConsecutives))

  tmp <- full_join(dats, tmp, by = "ID")

  return(tmp)
}

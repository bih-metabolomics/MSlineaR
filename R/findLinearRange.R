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
findLinearRange <- function(dat, x="DilutionPoint", y = "IntensityNorm", modelObject, res = 2, minConsecutives){
  #browser()

  setorder(dat,DilutionPoint)
  dat <- dat[color %in% "black"]
  modelObject <- unique(dat[[modelObject]])
  modelObject <- unlist(modelObject, recursive = F)
  int50 <- DescTools::Closest(x = dat[[y]] ,a = (max(modelObject$fit) -min(modelObject$fit))/2 + min(modelObject$fit), which = TRUE, na.rm = T)

  if(int50 == length(dat[[x]])) int50 <- length(dat[[x]]) -1
  if(int50 == 1) int50 = 2
  dat$color[int50] <-  "green"
  dat$pch[int50] <- 19



  #create linear regression line going through int50
  we <- rep(1, length(dat[[x]]))
  we[(int50 - 1) : (int50 + 1)] <- 1000


  linearRange <- lm(modelObject$fit ~ dat[[x]], weights = we)
  ablineIntensity <- fitted(linearRange)

  #confint <- linearRange$coefficients[1]*res
  #confint <- max(dat[[y]])/100*res
  confi <- abs(ablineIntensity/100)*res



  #if(all(ablineIntensity/100*(100 + res) >= ablineIntensity)){
    limit1 <- ablineIntensity + confi
    limit2 <- ablineIntensity - confi

    if(limit1[1] < limit2[1]){
      limitdown <- limit1
      limitup <- limit2
    } else{
      limitdown <- limit2
      limitup <- limit1
    }
    #limitup <- (linearRange$coefficients[1] + confint) + linearRange$coefficients[2]*dat[[x]]
    #limitdown <- (linearRange$coefficients[1] - confint) + linearRange$coefficients[2]*dat[[x]]
  #} else{

    #limitup <- ablineIntensity - confi
    #limitdown <- ablineIntensity + confi
    #limitup <- (linearRange$coefficients[1] - confint) + linearRange$coefficients[2]*dat[[x]]
    #limitdown <- (linearRange$coefficients[1] + confint) + linearRange$coefficients[2]*dat[[x]]

  #}

  consNDX <- rle(data.table::between(x = dat[[y]], lower = limitdown, upper = limitup, NAbounds = NA, check = T))




#rm(linearRange)
  #ndx <- which(abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))) < res)
  #consNDX <- rle(round(abs((dat[[y]] - ablineIntensity) /max(abs(dat[[y]] - ablineIntensity))), digits = 2) <= res)


  consNDX$position <- cumsum(consNDX$length)

  if(any(consNDX$length[which(consNDX$values %in% TRUE)]>= 3) & any(consNDX$position[consNDX$values %in% TRUE] >= int50)){

  TRUEpos <- Position(function(fi) fi >= int50, consNDX$position[consNDX$values %in% TRUE], right = TRUE)
  maxTrueRange <- (consNDX$position[consNDX$values %in% TRUE
                                   ][TRUEpos] - consNDX$length[consNDX$values %in% TRUE][TRUEpos] +1) : consNDX$position[consNDX$values %in% TRUE][TRUEpos]
maxTrueRange <- maxTrueRange[maxTrueRange!=0]
 # if(length(consNDX$length[which(consNDX$values %in% TRUE)]))



    tmp <- tibble(
      IDintern = dat$IDintern,
      groupIndices = dat$groupIndices,
      linearRangeStart = dat$DilutionPoint[maxTrueRange[1]],
      linearRangeEnd = dat$DilutionPoint[tail(maxTrueRange,1)],
      IslinearRange = dat$DilutionPoint >= linearRangeStart & dat$DilutionPoint <= linearRangeEnd,
      IsPositivAssociated = (dat[[y]] - lag(dat[[y]])) > 0,
      linearRange = length(maxTrueRange),
      enoughPointsWithinLinearRange = linearRange >= minConsecutives,
      color = dat$color,
      pch = dat$pch,
      modelFit = modelObject$fit,
      ablineFit = ablineIntensity,
      ablineLimit1 = limitup,
      ablineLimit2 = limitdown,
      ablineRes = (dat[[y]] - ablineIntensity) /dat[[y]],
      Comment = "linearRange found"

    )

    if(any(tmp$IsPositivAssociated[tmp$linearRangeStart : tmp$linearRangeEnd] %in% FALSE)){
      tmp$IslinearRange[tmp$IsPositivAssociated %in% FALSE] = FALSE
      countT <- rle(tmp$IslinearRange)
      countT$position <- cumsum(countT$length)

      if(any(countT$length[which(countT$values %in% TRUE)]>= 3) & any(countT$position[countT$values %in% TRUE] >= int50)){

        TRUEpos <- Position(function(fi) fi >= int50, countT$position[countT$values %in% TRUE], right = TRUE)
        maxTrueRange <- (countT$position[countT$values %in% TRUE
        ][TRUEpos] - countT$length[countT$values %in% TRUE][TRUEpos] +1) : countT$position[countT$values %in% TRUE][TRUEpos]
        maxTrueRange <- maxTrueRange[maxTrueRange!=0]

        tmp$linearRangeStart = dat$DilutionPoint[maxTrueRange[1]]
        tmp$linearRangeEnd = dat$DilutionPoint[tail(maxTrueRange,1)]
        tmp$IslinearRange = dat$DilutionPoint >= tmp$linearRangeStart & dat$DilutionPoint <= tmp$linearRangeEnd
        tmp$linearRange = length(maxTrueRange)
        tmp$enoughPointsWithinLinearRange = tmp$linearRange >= minConsecutives
      } else{

      tmp$IslinearRange = FALSE
      tmp$Comment = "linearRange found, but not positive associated"
      tmp$enoughPointsWithinLinearRange = FALSE
    }

  }} else{
    tmp <- tibble(
      IDintern = dat$IDintern,
      groupIndices = dat$groupIndices,
      linearRangeStart = NA,
      linearRangeEnd = NA,
      IslinearRange = FALSE,
      IsPositivAssociated = (dat[[y]] - lag(dat[[y]])) > 0,
      linearRange = NA,
      enoughPointsWithinLinearRange = NA,
      color = dat$color,
      pch = dat$pch,
      modelFit = modelObject$fit,
      ablineFit = ablineIntensity,
      ablineLimit1 = limitup,
      ablineLimit2 = limitdown,
      ablineRes = (dat[[y]] - ablineIntensity) /dat[[y]],
      Comment = "no linearRange found"

    )
  }
  if(tmp$IsPositivAssociated[2] %in% TRUE) tmp$IsPositivAssociated[1] <- TRUE

  return(tmp)
}



#--------------------------------------------------------------




# function func_lm performs linear regression for 3 consecutive points, starting at point x
# input : data table with column names [Intensity, Concentration, ID]
# output: tibble with informations about ID, analyzed point x, n distance and degree, and adjusted R2
#' Title
#'
#' @param tbl
#' @param x
#' @param explanVar
#'
#' @return
#' @export
#'
#' @examples
func_lm <- function(tbl, x, explanVar = 'Concentration'){

  tbl_linear_Range <- list()
  lm_sum <- lm(Intensity ~ get(explanVar), data = tbl[x : (x + 2),]) %>%
    summary()

  slope <- coef(lm_sum)[2]
  alpha.rad <- atan(slope)

  alpha.deg <- 180*alpha.rad/pi


  tibble(
    metabolite = tbl$ID[x],
    seqStart = tbl$DilutionPoint[x],
    seqEnd = tbl$DilutionPoint[x + 2],
    slope = lm_sum$coefficient[2],
    adjr2 = lm_sum$adj.r.squared,
    degree = alpha.deg#,
    #indice = tbl$indice[x]
  )
}

# function func_findLinearRange compares adjusted R2 and slope in degree to defined borders
# input : data table from func_lm with column names [metabolite, adjR2, degree]; values for adjusted R2 and degree; Number of points necessary to consider a linear range
# output: tibble with informations about ID, start of linear range, end of linear range, Comment
func_findLinearRange <- function(tbl, adjR2 , degr, points){

  adjr2.compare <- (round(tbl$adjr2,2) >= adjR2) & (tbl$degree >= degr)

  if (!any(adjr2.compare)) { # adjusted R2 and /or slope are to small
    tibble(
      ID = unique(tbl$metabolite),
      start = NA,
      end = NA,
      lengthRange = NA,
      Comment = "no linear range"
    )

  }else{ # minimum of 3 consecutive points with sufficient adjusted R2 and slope


    cons <- rle(adjr2.compare) # Compute the lengths of consecutive Trues and Falses
    indic <- which(cons$values == TRUE & cons$lengths >= (points - 2) )

    if (length(indic) > 0 ){ # length of consecutive Trues are sufficent
      cons.lengths.cumsum.start = cumsum(cons$lengths) - cons$lengths + 1
      cons.lengths.cumsum.end = cumsum(cons$lengths)
      starts = tbl$seqStart[cons.lengths.cumsum.start[indic]]
      ends = tbl$seqEnd[cons.lengths.cumsum.end[indic]]
      maxlin  <- tbl[which(tbl[ , "degree"] == max(tbl[cons.lengths.cumsum.start[indic]:cons.lengths.cumsum.end[indic], "degree"])), ]
      #newindex = ifelse(indic > 1, indic - 1, 0)
      #starts = cons.lengths.cumsum[newindex] + 1
      #if (0 %in% newindex) starts = c(1,starts)

      tibble(
        ID = unique(tbl$metabolite),
        start = starts,
        end = ends,
        lengthRange = ends - starts + 1,
        maxLinearRange = paste(maxlin[, c("seqStart", "seqEnd")], collapse = " : "),
        maxadjR2 = max(tbl[cons.lengths.cumsum.start[indic]:cons.lengths.cumsum.end[indic], "adjr2"]),
        Comment = ""
      )
    } else { # length of consecutive Trues are to small

      tibble(
        ID = unique(tbl$metabolite),
        start = NA,
        end = NA,
        lengthRange = NA,
        maxLinearRange = NA,
        maxadjR2 = NA,
        Comment = "no linear range"
      )
    }


  }
  #})
}

# function func_slider_LinearRange slides through data table of metabolites
# input: data table with column names [ID, Intensity, Concentration]; Number of points necessary to consider a linear range, boundaries for adjusted R2 value and slope in degree; cooks distance
# output:

func_slider_LinearRange <- function(dat, metabolite, points = 5, adjR2 = 0.85, degr = 10, residual = 2){ #, cooksdis = 4

  print(metabolite)
  tbl_linear_Range <- list()
  linRange <- data.frame()
  data <- dat %>% filter(ID == metabolite) %>% arrange(Concentration) %>% mutate(DilutionPoint = 1 : nrow(.), Comment = "")
  dataRaw <- data

  #outlier
  data <- data[!is.na(data$Intensity), ]
  model <- lm(Intensity ~ Concentration, data)
  data$residuals <- residuals(model)/sd(residuals(model))
  outlier <- which(abs(data$residuals) > residual)

  dataRaw <- join(type = "right", data, dataRaw)
  dataRaw$Comment[abs(dataRaw$residuals) > residual] = "Outlier"
  dataRaw$Comment[is.na(dataRaw$Intensity)] = "Missing"
  tbl_linear_Range[[metabolite]][["rawData"]] <- dataRaw

  if (any(outlier) | any(is.na(dataRaw$Intensity))) {
    # outlier
    nrOutlier <- length(outlier)
    linRange <- tibble()
    if (sum(outlier) > 0) data <- dataRaw[-outlier, ]
    # NA
    nrNA <- sum(is.na(dataRaw$Intensity))
    NApos <- which(is.na(dataRaw$Intensity))
    linRange <- tibble()
    data <- data[!is.na(data$Intensity), ]

  }



  tbl_lm <- slide(1:(nrow(data) - 2), ~func_lm(tbl = data, .x )) %>% ldply

  tbl_linear_Range[[metabolite]][["allLinearRanges"]] <- tbl_lm


  linRange <- func_findLinearRange(tbl = tbl_lm, adjR2, degr, points )
  linRange <- linRange %>% filter(lengthRange == max(lengthRange))

  if (any(linRange$Comment == "")) { # found one linear range, calculate slope and R2 for the whole range
    sum.all <- lm(Intensity ~ Concentration, data = dataRaw[linRange$start:linRange$end,]) %>%
      summary()
    alpha.rad <- atan(
      (dataRaw$Intensity[linRange$end] - dataRaw$Intensity[linRange$start])/
        (dataRaw$Concentration[linRange$end] - dataRaw$Concentration[linRange$start]))
    alpha.deg <- 180*alpha.rad/pi


    linRange$IntensityStart <- dataRaw$Intensity[linRange$start]
    linRange$IntensityEnd <- dataRaw$Intensity[linRange$end]
    linRange$ConcenstrationStart <- dataRaw$Concentration[linRange$start]
    linRange$ConcenstrationEnd <- dataRaw$Concentration[linRange$end]
    linRange$LinearRangeStart <- linRange$start
    linRange$LinearRangeEnd <- linRange$end
    linRange$adjr2 <- sum.all$adj.r.squared
    linRange$degree <- alpha.deg
    linRange$maxLinearRange <- linRange$maxLinearRange
    linRange$maxadjR2 <- linRange$maxadjR2
    linRange$Comment <- paste0(sum(dataRaw$Comment == "Outlier"), " outliers and ", sum(dataRaw$Comment == "Missing"), " Missings excluded")
    linRange$PositionOfMissings <- ifelse(any(is.na(dataRaw$Intensity)),paste(dataRaw$DilutionPoint[dataRaw$Comment == "Missing"], collapse = "; "), NA)
    linRange$PositionOfOutliers <- ifelse(any(dataRaw$Comment %in% "Outlier"),paste(dataRaw$DilutionPoint[dataRaw$Comment == "Outlier"], collapse = "; "), NA)
    linRange$OutlierInLinearRange <- ifelse(any(dataRaw$Comment %in% "Outlier"),any(between(outlier,linRange$start, linRange$end )),NA)
    linRange$MissingInLinearRange <- ifelse(any(is.na(dataRaw$Intensity)),any(between(NApos,linRange$start, linRange$end )), NA)
    linRange$type <- case_when(
      identical(sort(unique(dataRaw$Comment)), c("","Outlier")) ~ "Outlier",
      identical(sort(unique(dataRaw$Comment)), c("","Missing")) ~ "Missing",
      identical(sort(unique(dataRaw$Comment)), c("","Missing","Outlier")) ~ c("Outlier ,Missing"),
      TRUE ~ "completeCase"
    )


  } else {# no linear range found
    linRange <- tibble(ID = metabolite)

    linRange$IntensityStart <- NA
    linRange$IntensityEnd <- NA
    linRange$ConcenstrationStart <- NA
    linRange$ConcenstrationEnd <- NA
    linRange$LinearRangeStart <- NA
    linRange$LinearRangeEnd <- NA
    linRange$lengthRange <- NA
    linRange$adjr2 <- NA
    maxLinearRange <- NA
    linRange$degree <- NA
    linRange$maxadjR2 <- NA
    linRange$maxLinearRange <- NA
    linRange$Comment <- "no linear range"
    linRange$PositionOfMissings <- ifelse(any(is.na(dataRaw$Intensity)),paste(dataRaw$DilutionPoint[dataRaw$Comment == "Missing"]), "")
    linRange$PositionOfOutliers <- ifelse(any(dataRaw$Comment %in% "Outlier"),paste(dataRaw$DilutionPoint[dataRaw$Comment == "Outlier"]), "")
    linRange$OutlierInLinearRange <- NA
    linRange$MissingInLinearRange <- NA
    linRange$type <- case_when(
      identical(sort(unique(dataRaw$Comment)), c("","Outlier")) ~ "Outlier",
      identical(sort(unique(dataRaw$Comment)), c("","Missing")) ~ "Missing",
      identical(sort(unique(dataRaw$Comment)), c("","Missing","Outlier")) ~ c("Outlier ,Missing"),
      TRUE ~ "completeCase"
    )

  }

  tbl_linear_Range[[metabolite]][["SummaryLinearRange"]] <- linRange %>% select(ID, IntensityStart, IntensityEnd, ConcenstrationStart, ConcenstrationEnd, LinearRangeStart, LinearRangeEnd, lengthRange, adjr2, degree, maxadjR2,maxLinearRange, Comment, PositionOfOutliers, OutlierInLinearRange, PositionOfMissings, MissingInLinearRange, type)

  tbl_linear_Range[[metabolite]]$rawData <-
    tbl_linear_Range[[metabolite]]$rawData %>%
    filter(between(
      DilutionPoint,
      tbl_linear_Range[[metabolite]][["SummaryLinearRange"]]$LinearRangeStart,
      tbl_linear_Range[[metabolite]][["SummaryLinearRange"]]$LinearRangeEnd)) %>%
    mutate(Comment = replace(Comment, row_number() == 1, "start"), Comment = replace(Comment, row_number() == n(), "end")) %>%

    mutate(Comment = ifelse(Comment == "", "linear", Comment)) %>%
    rbind(tbl_linear_Range[[metabolite]]$rawData %>%
            filter( DilutionPoint < tbl_linear_Range[[metabolite]][["SummaryLinearRange"]]$LinearRangeStart |
                      DilutionPoint > tbl_linear_Range[[metabolite]][["SummaryLinearRange"]]$LinearRangeEnd)) %>%
    arrange(DilutionPoint)

  return(tbl_linear_Range)
}








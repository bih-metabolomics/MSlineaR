
#' Title
#'
#' @param LR_object
#' @param printPDF
#' @param groupIndices
#' @param statusLinear
#' @param ID
#' @param outputfileName
#' @param pdfwidth
#' @param pdfheight
#'
#' @return
#' @export
#'
#' @examples
plotFDS <- function(LR_object, printPDF = TRUE, groupIndices = "all", statusLinear = c(TRUE, FALSE), ID = "all", outputfileName = c("Calibrationplot"), pdfwidth = 9, pdfheight = 60){

  plotDataGroup <- data.table(LR_object[["summaryFDS"]])
  plotDataGroup$enoughPointsWithinLR[is.na(plotDataGroup$enoughPointsWithinLR)] <- FALSE
  plotDataGroup <- plotDataGroup[enoughPointsWithinLR %in% statusLinear]
  if(groupIndices != "all" & groupIndices != "") plotDataGroup <- plotDataGroup[groupIndices]
  if(ID != "all" & ID != "") plotDataGroup <- plotDataGroup[get(LR_object$Parameters$COLNAMES[["ID"]]) %in% ID]

  plotDataFeature <- LR_object[["summaryFFDS"]]
  plotDataFeature <- plotDataFeature[groupIndices %in% plotDataGroup$groupIndices, ]

  plotDataFeature <- plotDataFeature[, ':=' (plotID = get(LR_object$Parameters$COLNAMES[["ID"]]),
                                             Yvalue = get(LR_object$Parameters$Y),
                                             Xvalue = get(LR_object$Parameters$X),
                                             BatchVar = get(LR_object$Parameters$COLNAMES[["REPLICATE"]]))]


  plotDataFeature <- plotDataFeature[plotDataGroup[,.(groupIndices, R2)], on = "groupIndices"]

  plotlinearData <- plotDataFeature |>
    ggplot(mapping = aes(x = Xvalue, y = Yvalue, col = color)) +
    geom_point(aes(color = color)) +
    scale_shape_identity() +
    scale_colour_identity() +
    geom_line(aes(x = Xvalue, y = abline), col = "orange") +
    geom_line(aes(x = Xvalue, y = ablineLimit1), linetype = "dotted", col = "red") +
    geom_line(aes(x = Xvalue, y = ablineLimit2), linetype = "dotted", col = "red") +
    geom_vline(data = plotDataGroup,aes( xintercept = log(plotDataGroup$LRStartX), color = "darkgrey"), linetype = "dotted") +
    geom_vline(data = plotDataGroup, aes( xintercept = log(plotDataGroup$LREndX), color = "darkgrey"), linetype = "dotted") +

    facet_grid(get(LR_object$Parameters$COLNAMES[["ID"]]) ~ get(LR_object$Parameters$COLNAMES[["REPLICATE"]]) ,  scales = "free") +
    geom_text(aes(x = -Inf, y = Inf, label = paste("R2 = ", round(R2,2)), group = Compound),
              size = 5,
              hjust = -0.5,
              vjust = 1.4,
              inherit.aes = FALSE) +
    theme_bw() +
    theme(panel.grid.minor=element_blank()) +
    theme(panel.grid.major=element_blank()) +
    theme(panel.background=element_blank()) +
    theme(axis.line=element_line()) +
    theme(legend.position="top")

  if(LR_object$Parameters$LOG_TRANSFORM %in% TRUE) plotlinearData <- plotlinearData + scale_y_continuous(name = "Area", labels = scales::trans_format("exp"))



  if(printPDF %in% TRUE){

  pdf(file = file.path(getwd(),paste0(outputfileName,".pdf")), width = pdfwidth, height = pdfheight)
  plot(plotlinearData)
  dev.off()
  } else {
    plotlinearData
  }
  return(plotlinearData)

}


#' Title
#'
#' @param LR_object
#' @param data_Sample
#' @param groupIndices
#' @param statusLinear
#' @param ID
#' @param printPDF
#' @param outputfileName
#' @param pdfwidth
#' @param pdfheight
#'
#' @return
#' @export
#'
#' @examples
plotSamples <- function(LR_object, data_Sample, groupIndices = "all", statusLinear = c(TRUE, FALSE), ID = "all", printPDF = TRUE, outputfileName = c("Sampleplot"), pdfwidth = 9, pdfheight = 60){


  plotCalibrant <- plotFDS(LR_object,
          groupIndices,
          statusLinear,
          ID,
          printPDF = FALSE)

  plotDataGroup <- data.table(LR_object[["summaryFDS"]])

  data_tbl_Sample <- getConc(dats = data_Sample, plotDataGroup)

  plotDataSample <- data_tbl_Sample[get(LR_object$Parameters$COLNAMES[["ID"]]) %in% plotCalibrant$data$plotID]
  plotDataGroup <- plotDataGroup[get(LR_object$Parameters$COLNAMES[["ID"]]) %in% plotCalibrant$data$plotID]


  if(LR_object$Parameters$LOG_TRANSFORM %in% TRUE){
    plotSamples <- plotCalibrant +
      geom_point(data = plotDataSample[Status_LR %in% TRUE], aes( x = log(ConcentrationLR), y = log(get(LR_object$Parameters$COLNAMES[["Y"]])), color = "purple"), inherit.aes = FALSE) +
      geom_point(data = plotDataSample[Status_LR %in% FALSE], aes( x = log(ConcentrationLR), y = log(get(LR_object$Parameters$COLNAMES[["Y"]])), color = "red"), inherit.aes = FALSE)

  } else{


    plotSamples <- plotCalibrant +
      geom_point(data = plotDataSample[Status_LR %in% TRUE], aes( x = ConcentrationLR, y = get(LR_object$Parameters$COLNAMES[["Y"]]), color = "purple"), inherit.aes = FALSE) +
      geom_point(data = plotDataSample[Status_LR %in% FALSE], aes( x = ConcentrationLR, y = get(LR_object$Parameters$COLNAMES[["Y"]]), color = "red"), inherit.aes = FALSE) +
      geom_vline(data = plotDataGroup,aes( xintercept = plotDataGroup$LRStartX, color = "darkgrey"), linetype = "dotted", inherit.aes = FALSE) +
      geom_vline(data = plotDataGroup, aes( xintercept = plotDataGroup$LREndX, color = "darkgrey"), linetype = "dotted", inherit.aes = FALSE)
    }




  if(printPDF %in% TRUE){

    pdf(file = file.path(getwd(),paste0(outputfileName,".pdf")), width = pdfwidth, height = pdfheight)
    plot(plotSamples)
    dev.off()
  } else {
    plotSamples
  }
  return(plotSamples)

}







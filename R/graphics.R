
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
plotFDS <- function(inputData_Series, inputData_BioSamples, inputData_QC, inputData_Blank,
                    nrRow = 10,
                    printPDF = TRUE, GroupIndices = "all",  Feature = "all", printR2 = TRUE,
                    outputfileName = c("Calibrationplot"),
                    columns = c(ID = COLNAMES[["ID"]],Batch = COLNAMES[["Batch"]], X = Xraw, Y = Yraw), ...){
# LR_object,statusLinear = c(TRUE, FALSE),
  data_Signals <- inputData_Series
  ID <- columns[["ID"]]
  X <- columns[["X"]]
  Y <- columns[["Y"]]
  data.table::setDT(data_Signals)
  data.table::setDT(inputData_BioSamples)

  if(any(GroupIndices != "all" & GroupIndices != "")){
    data_Signals <- data_Signals[groupIndices %in% GroupIndices]
    inputData_BioSamples <- inputData_BioSamples[groupIndices %in% GroupIndices]
    inputData_QC <- inputData_QC[groupIndices %in% GroupIndices]
    }
  if(any(Feature != "all" & Feature != "")){
    data_Signals <- data_Signals[get(ID) %in% Feature]
    inputData_BioSamples <- inputData_BioSamples[get(ID) %in% Feature]
    inputData_QC <- inputData_QC[get(ID) %in% Feature]
    }


  data.table::setorderv(data_Signals, ID)
  data.table::setorderv(inputData_BioSamples, ID)

  nCol = data.table::uniqueN(data_Signals[,get(columns[["Batch"]])])
  npage = ceiling(data.table::uniqueN(data_Signals[,get(ID)])/nrRow)

  if(printPDF %in% TRUE){

    #plotObj <- vector("list", npage)

    pdf(file = file.path(getwd(),paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 9, height = 15)}
  for(page in 1:npage) {


  plotlinearData <- data_Signals |>
    ggplot2:: ggplot(mapping = aes(x = DilutionPoint, y = get(Y))) +
    geom_point(colour = "black")


  if("OutlierFOD" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      geom_point(data = subset(data_Signals, OutlierFOD %in% TRUE), aes(x = DilutionPoint, y = get(Y)), colour = "red")
  }

  if("trim" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      geom_point(data = subset(data_Signals, trim %in% TRUE), aes(x = DilutionPoint, y = get(Y)), colour = "grey")
  }

  if("OutlierSOD" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      geom_point(data = subset(data_Signals, OutlierSOD %in% TRUE), aes(x = DilutionPoint, y = get(Y)), colour = "red")
  }

  if("trimPos" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      geom_point(data = subset(data_Signals, trimPos %in% TRUE), aes(x = DilutionPoint, y = get(Y)), colour = "grey")
  }

  if("IsLinear" %in% colnames(data_Signals)){

      Xinterstart = data_Signals$DilutionPoint[data_Signals$LRStart]
      Xinterend = data_Signals$DilutionPoint[data_Signals$LREnd]


    if(TRANSFORM_Y !="" & !is.na(TRANSFORM_Y)){
      Yinterstart = get(TRANSFORM_Y)(data_Signals$LRStartY)
      Yinterend = get(TRANSFORM_Y)(data_Signals$LREndY)
    } else{
      Yinterstart =  data_Signals$LRStartY
      Yinterend = data_Signals$LREndY
    }


    plotlinearData <-  plotlinearData +
      geom_point(data = subset(data_Signals, IsLinear %in% TRUE), aes(x = DilutionPoint, y = get(Y)), colour = "seagreen") +
      geom_line(data = data_Signals,aes(x = DilutionPoint, y = abline), col = "orange") +
      geom_vline(data = data_Signals,aes( xintercept = Xinterstart), col = "darkgrey", linetype = "dotted") +
      geom_vline(data = data_Signals,aes( xintercept = Xinterend), col = "darkgrey", linetype = "dotted") +
      geom_hline(data = data_Signals,aes( yintercept = Yinterstart), col = "darkgrey", linetype = "dotted") +
      geom_hline(data = data_Signals,aes( yintercept = Yinterend), col = "darkgrey", linetype = "dotted")

  }




  if(!is.null(inputData_BioSamples )){
    plotlinearData <-  plotlinearData +
      scale_x_continuous(expand = c(0.1,0.2,0.1,0.1),breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
      geom_point(data = inputData_BioSamples[Status_LR %in% "TRUE"], aes( x =  min(data_Signals[[X]], na.rm = T) -3, y = get(Y), alpha = 0.5), shape = 21,col = "purple", size = 2) +
      geom_point(data = inputData_BioSamples[Status_LR %in% "FALSE" ], aes(x = min(data_Signals[[X]], na.rm = T) -3, y = get(Y), alpha = 0.5), shape = 21, col = "red", size = 2) +
      geom_vline(aes( xintercept = 0, color = "darkgrey"), linetype = "solid", col = "black") +
      geom_text(aes(x = -3, y = Inf, label = "Sample"), size = 3,vjust = 2)
  }

  if(!is.null(inputData_QC )){
    plotlinearData <-  plotlinearData +
      scale_x_continuous(expand = c(0.1,0.2,0.1,0.1),breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
      geom_point(data = inputData_QC[Status_LR %in% "TRUE"], aes( x = min(data_Signals[[X]], na.rm = T) -2, y = get(Y), alpha = 0.5), shape = 23, col = "purple", size = 2) +
      geom_point(data = inputData_QC[Status_LR %in% "FALSE" ], aes(x = min(data_Signals[[X]], na.rm = T) -2, y = get(Y), alpha = 0.5),shape = 23, col = "red", size = 2)
      #geom_vline(aes( xintercept = min(data_Signals[[X]], na.rm = T) -3, color = "darkgrey"), linetype = "solid", col = "black") +
      #geom_text(aes(x = min(data_Signals[[X]], na.rm = T) -5, y = Inf, label = "Sample"), size = 3,vjust = 2)
  }


  plotlinearData <- plotlinearData +

    #ylim(c(min(data_Signals[[Y]]-100), max(data_Signals[[Y]]) + 10000)) +
    #scale_x_continuous(name = "Concentration",breaks= unique(data_Signals[[X]]), labels = scales::trans_format(get(inverse_x))) +
    scale_y_continuous(name = "Area", labels = scales::trans_format(get(inverse_y)), limits = c(layer_scales(plotlinearData)$y$get_limits()[1], layer_scales(plotlinearData)$y$get_limits()[2] + 1 ))


  if(length(GroupIndices > 1) | GroupIndices %in% "all" | length(Feature > 1) | Feature %in% "all" ){
    plotlinearData <-  plotlinearData +
      ggforce::facet_grid_paginate(get(ID) ~ get(columns[["Batch"]]) ,  scales = "free", ncol = nCol,nrow = nrRow, page = page )
    }


  if(printR2 %in% TRUE) {plotlinearData <- plotlinearData +
    geom_text(data = subset(data_Signals, !is.na(R2)),
              aes(x = 1, y = Inf, label = paste(Series,": R2 = ", round(R2,2)) ,  group = get(ID)),
              size = 3,
              hjust = -0.1,
              vjust = 2,
              inherit.aes = FALSE)


  }






  plotlinearData <-  plotlinearData +
    theme_bw() +
    theme(panel.grid.minor=element_blank()) +
    theme(panel.grid.major=element_blank()) +
    theme(panel.background=element_blank()) +
    theme(axis.line=element_line()) +
    theme(axis.text.x = element_text(angle=90)) +
    theme(legend.position="none")


      plot(plotlinearData)

    }
  if(printPDF %in% TRUE){dev.off()}
  }



#   pdf(file = file.path(getwd(),paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 9, height = 25)
#    plot(plotlinearData)
#   dev.off()
#   } else {
#     plot(plotlinearData)
#   }
#   return(plotlinearData)
#
# }


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







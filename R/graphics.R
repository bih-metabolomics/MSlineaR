

#' Title
#'
#' @param inputData_Series
#' @param inputData_BioSamples
#' @param inputData_QC
#' @param inputData_QC_ref
#' @param inputData_Blank
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
combineData <- function(inputData_Series, inputData_BioSamples, inputData_QC, inputData_QC_ref, inputData_Blank
                        ){


  data_Signals <- inputData_Series

  if(!is.null(inputData_BioSamples)) data_Signals <- dplyr::full_join(data_Signals, inputData_BioSamples)
  if(!is.null(inputData_QC)) data_Signals <- dplyr::full_join(data_Signals, inputData_QC)
  if(!is.null(inputData_QC_ref)) data_Signals <- dplyr::full_join(data_Signals, inputData_QC_ref)
  if(!is.null(inputData_Blank)) data_Signals <- dplyr::full_join(data_Signals, inputData_Blank)




  return(data_Signals)

}




#' Title
#'
#' @param printPDF
#' @param inputData_Series
#' @param inputData_BioSamples
#' @param inputData_QC
#' @param inputData_QC_ref
#' @param inputData_Blank
#' @param nrRow
#' @param nrFeature
#' @param GroupIndices
#' @param Feature
#' @param printR2
#' @param ...
#' @param outputfileName
#'
#' @return
#' @export
#'
#' @examples
plot_FDS <- function(inputData_Series, inputData_BioSamples, inputData_QC, inputData_QC_ref = NULL, inputData_Blank = NULL,
                    nrRow = 10, nrFeature = 50,
                    printPDF = TRUE, GroupIndices = "all",  Feature = "all", printR2 = TRUE,
                    outputfileName = c("Calibrationplot"), TRANSFORM_Y, inverse_y,
                    COLNAMES, X, Y, Series ){

  assertthat::not_empty(inputData_Series)
# LR_object,statusLinear = c(TRUE, FALSE),

  ID <-  COLNAMES[["ID"]]
  Batch <-  COLNAMES[["Batch"]]
  X <- X
  Y <- Y

  data.table::setDT(inputData_Series)
  data.table::setDT(inputData_BioSamples)
  data.table::setDT(inputData_QC)
  data.table::setDT(inputData_QC_ref)
  data.table::setDT(inputData_Blank)



  data_Signals <- combineData(inputData_Series, inputData_BioSamples, inputData_QC, inputData_QC_ref, inputData_Blank)

  data.table::setorderv(data_Signals, ID)

  if(any(GroupIndices != "all" & GroupIndices != "")) data_Signals <- data_Signals[groupIndices %in% GroupIndices]
  if(any(Feature != "all" & Feature != "")) data_Signals <- data_Signals[get(ID) %in% Feature]
  if(Feature %in% "all" & GroupIndices %in% "all" & data.table::uniqueN(data_Signals[[ID]]) > nrFeature){
    randomIDs <- sample(unique(data_Signals[[ID]]), nrFeature, replace = F)
    data_Signals <- data_Signals[get(ID) %in% randomIDs]
  }




  nCol = data.table::uniqueN(data_Signals[[Batch]])
  npage = ceiling(as.numeric(data.table::uniqueN(data_Signals[[ID]])/nrRow))

  if(printPDF %in% TRUE){

    #plotObj <- vector("list", npage)

    pdf(file = file.path(getwd(),paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 9, height = 15)}
  for(page in 1:npage) {


  plotlinearData <-
    ggplot2:: ggplot(data = data_Signals, mapping = ggplot2::aes(x = DilutionPoint, y = get(Y)), shape = Sample.Type) +
    ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type])),colour = "black")


  if("OutlierFOD" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & OutlierFOD %in% TRUE), ggplot2::aes(x = DilutionPoint, y = get(Y)), colour = "red")
  }

  if("trim" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & trim %in% TRUE), ggplot2::aes(x = DilutionPoint, y = get(Y)), colour = "grey")
  }

  if("OutlierSOD" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & OutlierSOD %in% TRUE), ggplot2::aes(x = DilutionPoint, y = get(Y)), colour = "red")
  }

  if("trimPos" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & trimPos %in% TRUE), ggplot2::aes(x = DilutionPoint, y = get(Y)), colour = "grey")
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
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & IsLinear %in% TRUE), ggplot2::aes(x = DilutionPoint, y = get(Y)), colour = "seagreen") +
      ggplot2::geom_line(data = data_Signals,ggplot2::aes(x = DilutionPoint, y = abline), col = "orange") +
      ggplot2::geom_vline(data = data_Signals,ggplot2::aes( xintercept = Xinterstart), col = "darkgrey", linetype = "dotted") +
      ggplot2::geom_vline(data = data_Signals,ggplot2::aes( xintercept = Xinterend), col = "darkgrey", linetype = "dotted") +
      ggplot2::geom_hline(data = data_Signals,ggplot2::aes( yintercept = Yinterstart), col = "darkgrey", linetype = "dotted") +
      ggplot2::geom_hline(data = data_Signals,ggplot2::aes( yintercept = Yinterend), col = "darkgrey", linetype = "dotted")

  }

  plotlinearData <-  plotlinearData +
    ggplot2::scale_x_continuous(breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) #+#scales::trans_format(get(inverse_x), format = number_format())) +

# if(!is.null(inputData_BioSamples) | !is.null(inputData_QC)| !is.null(inputData_QCref) | !is.null(inputData_Blank)){
#   nrQC <- sum(!is.null(inputData_BioSamples),!is.null(inputData_QC),!is.null(inputData_QCref), !is.null(inputData_Blank))
#
# }


  if(!is.null(inputData_BioSamples )){
    plotlinearData <-  plotlinearData +
      ggplot2::scale_x_continuous(limits = c(-4, NA) ,breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_BioSamples[, Sample.Type]) ), ggplot2::aes( x = -1, y = get(Y),  shape = Sample.Type, color = Status_LR), size = 2) +#, shape = 1, col = "purple"
      ggplot2::geom_vline(ggplot2::aes( xintercept = 0, color = "darkgrey"), linetype = "solid", col = "black") +
      ggplot2::geom_text(ggplot2::aes(x = -2.5, y = Inf, label = "QC & Samples"), size = 3,vjust = 2)
    legend_order <- c(unique(inputData_BioSamples$Sample.Type))
  }

  if(!is.null(inputData_QC )){
    plotlinearData <-  plotlinearData +
      ggplot2::scale_x_continuous(limits = c(-4, NA) ,breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_QC[, Sample.Type])), ggplot2::aes( x = -2, y = get(Y),  shape = Sample.Type, color = Status_LR), size = 2) #+, shape = 5
      #geom_vline(ggplot2::aes( xintercept = min(data_Signals[[X]], na.rm = T) -3, color = "darkgrey"), linetype = "solid", col = "black") +
      #geom_text(ggplot2::aes(x = min(data_Signals[[X]], na.rm = T) -5, y = Inf, label = "Sample"), size = 3,vjust = 2)
    legend_order <- c(legend_order, unique(inputData_QC$Sample.Type))
  }

  if(!is.null(inputData_QC_ref )){
    plotlinearData <-  plotlinearData +
      ggplot2::scale_x_continuous(limits = c(-4, NA) ,breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_QC_ref[, Sample.Type])), ggplot2::aes( x = -3, y = get(Y),  shape = Sample.Type, color = Status_LR), size = 2) #+ shape = 2
    #geom_vline(ggplot2::aes( xintercept = min(data_Signals[[X]], na.rm = T) -3, color = "darkgrey"), linetype = "solid", col = "black") +
    #geom_text(ggplot2::aes(x = min(data_Signals[[X]], na.rm = T) -5, y = Inf, label = "Sample"), size = 3,vjust = 2)
    legend_order <- c(legend_order, unique(inputData_QC_ref$Sample.Type))

  }

  if(!is.null(inputData_Blank )){
    plotlinearData <-  plotlinearData +
      ggplot2::scale_x_continuous(limits = c(-4, NA) ,breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Blank[, Sample.Type])), ggplot2::aes( x = -4, y = get(Y),  shape = Sample.Type, color = Status_LR), size = 2)# + shape = 0
    #geom_vline(ggplot2::aes( xintercept = min(data_Signals[[X]], na.rm = T) -3, color = "darkgrey"), linetype = "solid", col = "black") +
    #geom_text(ggplot2::aes(x = min(data_Signals[[X]], na.rm = T) -5, y = Inf, label = "Sample"), size = 3,vjust = 2)
    legend_order <- c(legend_order, unique(inputData_Blank$Sample.Type))

  }



  plotlinearData <- plotlinearData +

    #ylim(c(min(data_Signals[[Y]]-100), max(data_Signals[[Y]]) + 10000)) +
    #scale_x_continuous(name = "Concentration",breaks= unique(data_Signals[[X]]), labels = scales::trans_format(get(inverse_x))) +
    ggplot2:: scale_y_continuous(name = "Area", labels = scales::trans_format(get(inverse_y)), limits = c(NA, ggplot2::layer_scales(plotlinearData)$y$get_limits()[2] + 1 ))


  if(length(GroupIndices > 1) | GroupIndices %in% "all" | length(Feature > 1) | Feature %in% "all" ){
    plotlinearData <-  plotlinearData +
      ggforce::facet_grid_paginate(ID ~ Batch ,  scales = "free", ncol = nCol,nrow = nrRow, page = page )
    }


  if(printR2 %in% TRUE) {plotlinearData <- plotlinearData +
    ggplot2::geom_text(data = subset(data_Signals, !is.na(R2)),
              ggplot2::aes(x = 0, y = Inf, label = paste(Series,": R2 = ", round(R2,2)) ,  group = get(ID)),
              size = 3,
              hjust = -0.1,
              vjust = 2,
              inherit.aes = FALSE)


  }






  plotlinearData <-  plotlinearData +
    ggplot2::scale_color_manual(name = "In linear Range:",
                       values = c("FALSE" = "red", "TRUE" = "purple")) +
    ggplot2::scale_shape_manual(values = c(0, 2, 5, 1)) +
    ggplot2::scale_shape_discrete(breaks=rev(legend_order)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
    ggplot2::theme(panel.background=ggplot2::element_blank()) +
    ggplot2::theme(axis.line=ggplot2::element_line()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90)) +
    ggplot2::theme(legend.position="top") +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
           shape = ggplot2::guide_legend(order = 2))


      plot(plotlinearData)

    }
  if(printPDF %in% TRUE){dev.off()}

  return(plotlinearData)
  }



#' Title
#'
#' @param printPDF
#' @param inputData_Series
#' @param GroupIndices
#' @param Feature
#' @param COLNAMES
#' @param X
#' @param Y
#' @param ...
#' @param outputfileName
#'
#' @return
#' @export
#'
#' @examples


plot_Barplot_Summary <- function(inputData_Series,
                                printPDF = TRUE, GroupIndices = "all",  Feature = "all",
                                outputfileName = c("Summary_Barplot"),
                                COLNAMES, X, Y , ...){

  assertthat::not_empty(inputData_Series)
  # LR_object,statusLinear = c(TRUE, FALSE),
  ID <- COLNAMES[["ID"]]
  X <- X
  Y <- Y

  data.table::setDT(inputData_Series)

  data_Signals_summary <- inputData_Series |>
    dplyr::group_by(DilutionPoint, Batch) |>
    dplyr::summarize(Missing = sum(is.na(Y), na.rm = T),
              OutlierFOD = sum(OutlierFOD, na.rm = T),
              Trim = sum(trim %in% TRUE, na.rm = T),
              OutlierSOD = sum(OutlierSOD, na.rm = T),
              LR_TRUE = sum(IsLinear %in% TRUE, na.rm = T),
              LR_FALSE = sum(IsLinear %in% FALSE, na.rm = T)
              ) |>
    tidyr::pivot_longer(names_to = "Type",values_to =  "count" ,cols = -c(1:2))



  data_Signals_summary$Type <- factor(data_Signals_summary$Type, levels = rev(c("LR_TRUE", "LR_FALSE", "Missing", "OutlierFOD", "OutlierSOD", "Trim")))
  data_Signals_summary <- data_Signals_summary |> dplyr::filter(count > 0)



  plot_Summary <-
    ggplot2::ggplot(data = data_Signals_summary, ggplot2::aes(x = DilutionPoint, y = count, label = count, fill = Type)) +
    ggplot2::geom_bar(stat="identity",
                 position="fill",
                 width = 0.5 ) +
    ggplot2::geom_text(size = 3, position = ggplot2::position_fill(vjust = 0.5)) +
    ggplot2::facet_grid(.~Batch, scales = "free_x", space = "free_x")

  plot_Summary <- plot_Summary +
    ggplot2::scale_x_continuous(breaks = data_Signals_summary$DilutionPoint) +
    ggplot2::scale_fill_discrete(name = "", ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_manual(values = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
    ggplot2::theme(panel.background=ggplot2::element_blank()) +
    ggplot2::theme(axis.line=ggplot2::element_line())

  if(printPDF %in% TRUE){

    pdf(file = file.path(getwd(),paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)}


  plot( plot_Summary)

  if(printPDF %in% TRUE){ dev.off()}

  return( plot_Summary)

}


#' Title
#'
#' @param inputData_Samples
#' @param printPDF
#' @param GroupIndices
#' @param Feature
#' @param outputfileName
#' @param ...
#' @param COLNAMES
#' @param X
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
plot_Barplot_Summary_Sample <- function(inputData_Samples,
                                        printPDF = TRUE, GroupIndices = "all",  Feature = "all",
                                        outputfileName = c("Summary_Barplot_Samples"),
                                        COLNAMES, X, Y , ...){


  assertthat::not_empty(inputData_Samples)
  # LR_object,statusLinear = c(TRUE, FALSE),
  ID <- COLNAMES[["ID"]]
  X <- X
  Y <- Y
  SAMPLE_ID <- COLNAMES[["Sample_ID"]]

  data.table::setDT(inputData_Samples)

  data_Signals_sample_summary <- inputData_Samples |>
    dplyr::group_by(Sample_ID = get(SAMPLE_ID), Batch, Sample.Type) |>
    dplyr::summarize(
      Missing = sum(is.na(Y), na.rm = T),
      LR_TRUE = sum(Status_LR %in% TRUE, na.rm = T),
      LR_FALSE = sum(Status_LR %in% FALSE, na.rm = T)
    ) |>
    tidyr::pivot_longer(names_to = "Type",values_to =  "count" ,cols = -c(1:3))



  data_Signals_sample_summary$Type <- factor(data_Signals_sample_summary$Type, levels = rev(c("LR_TRUE", "LR_FALSE", "Missing", "OutlierFOD", "OutlierSOD", "Trim")))
  data_Signals_sample_summary <- data_Signals_sample_summary |> dplyr::filter(count > 0)
  setorder(data_Signals_sample_summary, Sample.Type)
  data_Signals_sample_summary$Sample_ID <- factor(data_Signals_sample_summary$Sample_ID, levels = unique(data_Signals_sample_summary$Sample_ID))


  plot_Summary_samples <-
    ggplot2::ggplot(data = data_Signals_sample_summary, ggplot2::aes(x = Sample_ID, y = count, label = count, fill = Type)) +
    ggplot2::geom_bar(stat="identity",
             position="fill",
             width = 0.5 ) +
    #geom_text(size = 3, position = position_fill(vjust = 0.5)) +
    ggplot2::facet_grid(.~Batch, scales = "free_x", space = "free_x")

  plot_Summary_samples <- plot_Summary_samples +
    #scale_x_continuous(breaks = data_Signals$DilutionPoint) +
    ggplot2::scale_fill_discrete(name = "", ) +
    ggplot2:: scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_manual(values = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
    ggplot2:: theme(panel.background=ggplot2::element_blank()) +
    ggplot2:: theme(axis.text.x = ggplot2::element_text(angle=90)) +
    ggplot2:: theme(axis.line=ggplot2::element_line())

  if(printPDF %in% TRUE){

    pdf(file = file.path(getwd(),paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)}


  plot( plot_Summary_samples)

  if(printPDF %in% TRUE){ dev.off()}

  return( plot_Summary_samples)



}





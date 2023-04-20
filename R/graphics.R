

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
combineData <- function(inputData_Series, inputData_BioSamples, inputData_QC#, inputData_QC_ref, inputData_Blank
                        ){


  data_Signals <- inputData_Series

  if(!is.null(inputData_BioSamples)) data_Signals <- dplyr::full_join(data_Signals, inputData_BioSamples,
                                                                      by = intersect(colnames(data_Signals), colnames(inputData_BioSamples)))
  if(!is.null(inputData_QC)) data_Signals <- dplyr::full_join(data_Signals, inputData_QC,
                                                              by = intersect(colnames(data_Signals), colnames(inputData_QC)))
  # if(!is.null(inputData_QC_ref)) data_Signals <- dplyr::full_join(data_Signals, inputData_QC_ref,
  #                                                                 by = intersect(colnames(data_Signals), colnames(inputData_QC_ref)))
  # if(!is.null(inputData_Blank)) data_Signals <- dplyr::full_join(data_Signals, inputData_Blank,
  #                                                                by = intersect(colnames(data_Signals), colnames(inputData_Blank)))




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
plot_FDS <- function(inputData_Series, inputData_BioSamples, inputData_QC,
                     nrFeature = 50,
                    printPDF = TRUE, GroupIndices = "all",  Feature = "all", printR2 = TRUE,
                    outputfileName = c("Calibrationplot"), TRANSFORM_Y, inverse_y,
                    COLNAMES, X, Y, Series, output_dir ){

  assertthat::not_empty(inputData_Series)
# LR_object,statusLinear = c(TRUE, FALSE),

  ID <-  COLNAMES[["ID"]]
  Col_Batch <-  COLNAMES[["Batch"]]
  Sample.Type <- COLNAMES[["Sample_type"]]
  indipendent <- X
  y <- Y

  data.table::setDT(inputData_Series)
  data.table::setDT(inputData_BioSamples)
  data.table::setDT(inputData_QC)




  data_Signal <- combineData(inputData_Series, inputData_BioSamples, inputData_QC)#, inputData_QC_ref, inputData_Blank)
  data_Signal$ID = data_Signal[[ID]]
  data_Signal$Batch = data_Signal[[Col_Batch]]

  QCs = data.frame(Sample.Type = unique(inputData_QC[[Sample.Type]]),
                        x = -c(3 : (length(unique(inputData_QC[[Sample.Type]])) +2))
  )

  data_Signal <- dplyr::full_join(data_Signal, QCs, by = "Sample.Type")

  data.table::setorderv(data_Signal, ID)

  if(any(GroupIndices != "all" & GroupIndices != "")) data_Signal <- data_Signal[groupIndices %in% GroupIndices]
  if(any(Feature != "all" & Feature != "")) data_Signal <- data_Signal[get(ID) %in% Feature]
  if(Feature %in% "all" & GroupIndices %in% "all" & data.table::uniqueN(data_Signal[[ID]]) > nrFeature){
    randomIDs <- sample(unique(data_Signal[[ID]]), nrFeature, replace = F)
    data_Signal <- data_Signal[get(ID) %in% randomIDs]
  }





  if(printPDF %in% TRUE){

    #plotObj <- vector("list", npage)
    if(data.table::uniqueN(data_Signal$Batch) <=2){
      pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 9, height = 15)
      nrRow = 10
    }else {
      pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)
      nrRow = 5
    }
  }

    nCol = data.table::uniqueN(data_Signal[[Col_Batch]])
    npage = ceiling(as.numeric(data.table::uniqueN(data_Signal[[ID]])/nrRow))


  for(page in 1:npage) {

    data_Signals <- data_Signal#[get(ID) %in% unique(data_Signal[[ID]])[(nrRow * (page -1) +1) : (nrRow * page)]]

  plotlinearData <-
    ggplot2:: ggplot(data = data_Signals, mapping = ggplot2::aes(x = get(indipendent), y = get(y)), shape = Sample.Type) +
    ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type])),colour = "black")


  if("OutlierFOD" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & OutlierFOD %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(y)), colour = "red")
  }

  if("trim" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & trim %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(y)), colour = "grey")
  }

  if("OutlierSOD" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & OutlierSOD %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(y)), colour = "red")
  }

  if("trimPos" %in% colnames(data_Signals)){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & trimPos %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(y)), colour = "grey")
  }

  if("IsLinear" %in% colnames(data_Signals)){

      Xinterstart = data_Signals[[indipendent]][data_Signals$LRStart +1]
      Xinterend = data_Signals[[indipendent]][data_Signals$LREnd +1]


    if(TRANSFORM_Y !="" & !is.na(TRANSFORM_Y)){
      Yinterstart = get(TRANSFORM_Y)(data_Signals$LRStartY)
      Yinterend = get(TRANSFORM_Y)(data_Signals$LREndY)
    } else{
      Yinterstart =  data_Signals$LRStartY
      Yinterend = data_Signals$LREndY
    }


    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & IsLinear %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(y)), colour = "seagreen") +
      ggplot2::geom_line(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & IsLinear %in% TRUE),ggplot2::aes(x = get(indipendent), y = abline), col = "orange", na.rm = TRUE) +
      ggplot2::geom_vline(data = data_Signals,ggplot2::aes( xintercept = Xinterstart), col = "darkgrey", linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_vline(data = data_Signals,ggplot2::aes( xintercept = Xinterend), col = "darkgrey", linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_hline(data = data_Signals,ggplot2::aes( yintercept = Yinterstart), col = "darkgrey", linetype = "dotted", na.rm = TRUE) +
      ggplot2::geom_hline(data = data_Signals,ggplot2::aes( yintercept = Yinterend), col = "darkgrey", linetype = "dotted", na.rm = TRUE)

  }

  plotlinearData <-  plotlinearData +
    ggplot2::scale_x_continuous(name = "Dilution", limits = c(-5, NA) ,breaks = data_Signals[[indipendent]],  labels = data_Signals$DilutionPoint +1) #+#scales::trans_format(get(inverse_x), format = number_format())) +

# if(!is.null(inputData_BioSamples) | !is.null(inputData_QC)| !is.null(inputData_QCref) | !is.null(inputData_Blank)){
#   nrQC <- sum(!is.null(inputData_BioSamples),!is.null(inputData_QC),!is.null(inputData_QCref), !is.null(inputData_Blank))
#
# }


  if(!is.null(inputData_BioSamples )){
    plotlinearData <-  plotlinearData +
      #ggplot2::scale_x_continuous(limits = c(-4, NA) ,breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_BioSamples[, Sample.Type]) ), ggplot2::aes( x = -2, y = get(y),  shape = Sample.Type, color = Status_LR), size = 2, na.rm = TRUE) +#, shape = 1, col = "purple"
      ggplot2::geom_vline(ggplot2::aes( xintercept = -1, color = "darkgrey"), linetype = "solid", col = "black", na.rm = TRUE)
    legend_order <- c(unique(inputData_BioSamples$Sample.Type))
  }

  if(!is.null(inputData_QC )){
    plotlinearData <-  plotlinearData +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% QCs$Sample.Type),
                          ggplot2::aes( x = x, y = get(y),  shape = Sample.Type,
                                        color = Status_LR), size = 2, na.rm = TRUE) #+, shape = 5
    legend_order <- c(legend_order, unique(inputData_QC$Sample.Type))
  }



  plotlinearData <- plotlinearData +
    ggplot2:: scale_y_continuous(name = "Area", labels = scales::trans_format(get(inverse_y)),
                                 limits = c(NA, ggplot2::layer_scales(plotlinearData)$y$get_limits()[2] + 1 ))






  if(length(GroupIndices > 1) | GroupIndices %in% "all" | length(Feature > 1) | Feature %in% "all" ){
    plotlinearData <-  plotlinearData +
      ggforce::facet_grid_paginate(ID ~ Batch ,
                                   scales = "free", ncol = nCol,nrow = nrRow, page = page )
  }


  if(printR2 %in% TRUE) {

    text_label <- data_Signals[, .(ID, Batch, R2)]
    text_label <- text_label |> dplyr::group_by(ID, Batch) |>
      dplyr::mutate(R2 = ifelse(any(!is.na(R2)), R2[!is.na(R2)] , NA))
    text_label <- unique(text_label)




    plotlinearData <- plotlinearData +
      ggplot2::geom_text(data = text_label,
                         #ggplot2::geom_text(data = text_label[1:20,],#subset(data_Signals, !is.na(R2)),
                         ggplot2::aes(x = 0, y = Inf, label = paste(Series,": R2 = ", round(R2,4))),
                         size = 3,
                         hjust = -0.1,
                         vjust = 2,
                         #inherit.aes = FALSE
      )


  }







  plotlinearData <-  plotlinearData +
    ggplot2::geom_text(ggplot2::aes(x = -3.5, y = Inf, label = "QC & Samples"), size = 3,vjust = 2, na.rm = TRUE)+
    ggplot2::scale_color_manual(name = "In linear Range:",
                       values = c("FALSE" = "red", "TRUE" = "purple")) +
    ggplot2::scale_shape_manual(values = c(0, 2, 5, 1), breaks=rev(legend_order)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
    ggplot2::theme(panel.background=ggplot2::element_blank()) +
    ggplot2::theme(axis.line=ggplot2::element_line()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90)) +
    ggplot2::theme(legend.position="top") +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
           shape = ggplot2::guide_legend(order = 2))


      print(plotlinearData)

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
                                COLNAMES, X, Y , output_dir){

  assertthat::not_empty(inputData_Series)
  # LR_object,statusLinear = c(TRUE, FALSE),
  ID <- COLNAMES[["ID"]]
  Col_Batch = COLNAMES[["Batch"]]
  x <- X
  y <- Y

  data.table::setDT(inputData_Series)

  data_Signals_summary <- inputData_Series |>
    dplyr::group_by(DilutionPoint, Batch = get(Col_Batch)) |>
    dplyr::summarize(Missing = sum(is.na(y), na.rm = T),
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
    ggplot2::geom_text(size = 3, position = ggplot2::position_fill(vjust = 0.5), na.rm = TRUE) +
    ggplot2::facet_grid(. ~ Batch, scales = "free_x", space = "free_x")

  plot_Summary <- plot_Summary +
    ggplot2::scale_x_continuous(breaks = data_Signals_summary$DilutionPoint) +
    #ggplot2::scale_fill_discrete(name = "", ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_manual(name = "", values = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
    ggplot2::theme(panel.background=ggplot2::element_blank()) +
    ggplot2::theme(axis.line=ggplot2::element_line())

  if(printPDF %in% TRUE){

    pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)}


  print( plot_Summary)

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
                                        COLNAMES, X, Y , output_dir){


  assertthat::not_empty(inputData_Samples)
  # LR_object,statusLinear = c(TRUE, FALSE),
  ID <- COLNAMES[["ID"]]
  Col_Batch <- COLNAMES[["Batch"]]
  x <- X
  y <- Y
  SAMPLE_ID <- COLNAMES[["Sample_ID"]]

  data.table::setDT(inputData_Samples)

  data_Signals_sample_summary <- inputData_Samples |>
    dplyr::group_by(Sample_ID = get(SAMPLE_ID), Batch = get(Col_Batch), Sample.Type) |>
    dplyr::summarize(
      Missing = sum(is.na(y), na.rm = T),
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
    ggplot2::facet_grid(. ~ Batch, scales = "free_x", space = "free_x")

  plot_Summary_samples <- plot_Summary_samples +
    #scale_x_continuous(breaks = data_Signals$DilutionPoint) +
    ggplot2:: scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_manual(name = "", values = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
    ggplot2:: theme(panel.background=ggplot2::element_blank()) +
    ggplot2:: theme(axis.text.x = ggplot2::element_text(angle=90)) +
    ggplot2:: theme(axis.line=ggplot2::element_line())

  if(printPDF %in% TRUE){

    pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)}


  plot( plot_Summary_samples)

  if(printPDF %in% TRUE){ dev.off()}

  return( plot_Summary_samples)



}







#' Title
#'
#' @param inputData_Series
#' @param inputData_BioSamples
#' @param inputData_QC
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
combineData <- function(inputData_Series, inputData_BioSamples, inputData_QC #inputData_Blank#, inputData_QC_ref, inputData_Blank
                        ){


  data_Signals <- inputData_Series

  if(!is.null(inputData_BioSamples)) data_Signals <- dplyr::full_join(data_Signals, inputData_BioSamples,
                                                                      by = intersect(colnames(data_Signals), colnames(inputData_BioSamples)))
  if(!is.null(inputData_QC)) data_Signals <- dplyr::full_join(data_Signals, inputData_QC,
                                                              by = intersect(colnames(data_Signals), colnames(inputData_QC)))
  # if(!is.null(inputData_QC_ref)) data_Signals <- dplyr::full_join(data_Signals, inputData_QC_ref,
  #                                                                 by = intersect(colnames(data_Signals), colnames(inputData_QC_ref)))
  #if(!is.null(inputData_Blank)) data_Signals <- dplyr::full_join(data_Signals, inputData_Blank,
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
plot_FDS <- function(inputData_Series, inputData_BioSamples, inputData_QC,#inputData_Blank,
                     nrFeature = 50,
                     signal_blank_ratio = 5,
                     printPDF = TRUE, GroupIndices = "all",  Feature = "all", printR2 = TRUE,
                     outputfileName = c("Calibrationplot"), TRANSFORM_Y = "log", inverse_y = "exp",
                     COLNAMES, Xcol, Ycol, Series = "QC dilution Curves", output_dir ){

  assertthat::not_empty(inputData_Series)
  # LR_object,statusLinear = c(TRUE, FALSE),

  ID <-  COLNAMES[["Feature_ID"]]
  Col_Batch <-  COLNAMES[["Batch"]]
  Sample.Type <- COLNAMES[["Sample_type"]]
  indipendent <- Xcol
  dependent <- Ycol
  ClassCol <- COLNAMES[["Class"]]

  data.table::setDT(inputData_Series)
  inputData_Series$Sample.Type <- inputData_Series[[Sample.Type]]

  if(!is.null(inputData_BioSamples)){

    inputData_BioSamples$Sample.Type <- inputData_BioSamples[[Sample.Type]]
    data.table::setDT(inputData_BioSamples)
    if(!is.null(TRANSFORM_Y)) {
      inputData_BioSamples$Sample_area <- get(TRANSFORM_Y)(inputData_BioSamples[[COLNAMES[["Y"]]]])
      data.table::setnames(inputData_BioSamples, old = "Sample_area",new = dependent)
      inputData_BioSamples <- inputData_BioSamples |> dplyr::select(which(!duplicated(names(inputData_BioSamples))))
    }

  }
  if(!is.null(inputData_QC)){
    inputData_QC$Sample.Type <- inputData_QC[[Sample.Type]]
    data.table::setDT(inputData_QC)
    if(!is.null(TRANSFORM_Y)) {
      inputData_QC$Sample_area <- get(TRANSFORM_Y)(inputData_QC[[COLNAMES[["Y"]]]])
      data.table::setnames(inputData_QC, old = "Sample_area",new = dependent)
      inputData_QC <- inputData_QC |> dplyr::select(which(!duplicated(names(inputData_QC))))

    }

  }
  #data.table::setDT(inputData_Blank)




  data_Signal <- combineData(inputData_Series, inputData_BioSamples, inputData_QC)# inputData_Blank)#, inputData_QC_ref, inputData_Blank)
  data_Signal$ID = data_Signal[[ID]]
  data_Signal$Batch = data_Signal[[Col_Batch]]
  data_Signal$Sample.Type = data_Signal[[Sample.Type]]


  if(!is.null(inputData_QC)){
    QCs = data.frame(Sample.Type = unique(inputData_QC[[Sample.Type]]),
                     x = -c(3 : (length(unique(inputData_QC[[Sample.Type]])) + 2)))





    data_Signal <- dplyr::full_join(QCs, data_Signal, by = "Sample.Type")

  } else {data_Signal$x = 0}

  data.table::setorderv(data_Signal, ID)
  data.table::setDT(data_Signal)

  if(any(GroupIndices != "all" & GroupIndices != "")) data_Signal <- data_Signal[groupIndices %in% GroupIndices]
  if(any(Feature != "all" & Feature != "")) data_Signal <- data_Signal[get(ID) %in% Feature]
  if(any(Feature %in% "all") & any(GroupIndices %in% "all") & data.table::uniqueN(data_Signal[[ID]]) > nrFeature){

    if(data.table::uniqueN(data_Signal[[ID]][data_Signal$Status_LR %in% TRUE]) > nrFeature/2) {
      randomIDsTRUE <- sample(unique(data_Signal[[ID]][data_Signal$Status_LR %in% TRUE]), nrFeature/2, replace = F)
      randomIDs2 <- sample(setdiff(unique(data_Signal[[ID]]), randomIDsTRUE), nrFeature/2, replace = F)
      randomIDs <- c(randomIDsTRUE, randomIDs2)
      data_Signal <- data_Signal[ID %in% randomIDs]
    } else if(data.table::uniqueN(data_Signal[[ID]][data_Signal$Status_LR %in% TRUE]) < nrFeature/2){
      randomIDsTRUE <- unique(data_Signal[[ID]][data_Signal$Status_LR %in% TRUE])
      randomIDs2 <- sample(setdiff(unique(data_Signal[[ID]]), randomIDsTRUE), nrFeature -length(randomIDsTRUE), replace = F)
      randomIDs <- c(randomIDsTRUE, randomIDs2)
      data_Signal <- data_Signal[ID %in% randomIDs]
    }
  }else if(any(Feature %in% "all") & any(GroupIndices %in% "all") & data.table::uniqueN(data_Signal[[ID]]) <= nrFeature){
    nrFeature = data.table::uniqueN(data_Signal[[ID]])

  }






    nrRow <- ifelse(data.table::uniqueN(data_Signal$Batch) <=2, 10, 5)
    nCol = data.table::uniqueN(data_Signal[[Col_Batch]])
    #npage = ceiling(as.numeric(data.table::uniqueN(data_Signal[[ID]])/nrRow))



    data_Signals <- data_Signal |> dplyr::arrange(groupIndices)#[groupIndices %in% 1, ]#[get(ID) %in% unique(data_Signal[[ID]])[(nrRow * (page -1) +1) : (nrRow * page)]]

    plotlinearData <-
      ggplot2:: ggplot(data = data_Signals, mapping = ggplot2::aes(x = get(indipendent), y = get(dependent)), shape = Sample.Type) +
      ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type])),colour = "black")

    if("signalBlankRatio" %in% colnames(data_Signals)){
      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) &
                                            signalBlankRatio %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(dependent)), colour = "grey")
    }
    if("OutlierFOD" %in% colnames(data_Signals)){
      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & OutlierFOD %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(dependent)), colour = "red")
    }

    if("trim" %in% colnames(data_Signals)){
      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & trim %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(dependent)), colour = "grey")
    }

    if("OutlierSOD" %in% colnames(data_Signals)){
      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & OutlierSOD %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(dependent)), colour = "red")
    }

    if("trimPos" %in% colnames(data_Signals)){
      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & trimPos %in% TRUE), ggplot2::aes(x = get(indipendent), y = get(dependent)), colour = "grey")
    }

    if("IsLinear" %in% colnames(data_Signals)){

      Xinterstart = data_Signals[[indipendent]][data_Signals$LRStart ]
      Xinterend = data_Signals[[indipendent]][data_Signals$LREnd ]


      if(!is.null(TRANSFORM_Y)){
        Yinterstart = get(TRANSFORM_Y)(data_Signals$LRStartY)
        Yinterend = get(TRANSFORM_Y)(data_Signals$LREndY)
      } else{
        Yinterstart =  data_Signals$LRStartY
        Yinterend = data_Signals$LREndY
      }


      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & IsLinear %in% TRUE),
                            ggplot2::aes(x = get(indipendent), y = get(dependent)), colour = "seagreen") +
        ggplot2::geom_line(data = subset(data_Signals, Sample.Type %in% unique(inputData_Series[, Sample.Type]) & IsLinear %in% TRUE),
                           ggplot2::aes(x = get(indipendent), y = abline), col = "orange", na.rm = TRUE) +
        ggplot2::geom_vline(data = data_Signals,ggplot2::aes( xintercept = Xinterstart), col = "darkgrey", linetype = "dotted", na.rm = TRUE) +
        ggplot2::geom_vline(data = data_Signals,ggplot2::aes( xintercept = Xinterend), col = "darkgrey", linetype = "dotted", na.rm = TRUE) +
        ggplot2::geom_hline(data = data_Signals,ggplot2::aes( yintercept = Yinterstart), col = "darkgrey", linetype = "dotted", na.rm = TRUE) +
        ggplot2::geom_hline(data = data_Signals,ggplot2::aes( yintercept = Yinterend), col = "darkgrey", linetype = "dotted", na.rm = TRUE)

    }

    plotlinearData <-  plotlinearData +
      ggplot2::scale_x_continuous(name = "Dilution", limits = c(min(data_Signals$x), NA) ,breaks = data_Signals[[indipendent]],  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
      ggplot2::geom_vline(ggplot2::aes( xintercept = -1, color = "darkgrey"), linetype = "solid", col = "black", na.rm = TRUE)

    # if(!is.null(inputData_BioSamples) | !is.null(inputData_QC)| !is.null(inputData_QCref) | !is.null(inputData_Blank)){
    #   nrQC <- sum(!is.null(inputData_BioSamples),!is.null(inputData_QC),!is.null(inputData_QCref), !is.null(inputData_Blank))
    #
    # }

    legend_order <- c()

    if(!is.null(inputData_BioSamples )){
      plotlinearData <-  plotlinearData +
        #ggplot2::scale_x_continuous(limits = c(-4, NA) ,breaks = data_Signals$DilutionPoint,  labels = data_Signals$DilutionPoint) +#scales::trans_format(get(inverse_x), format = number_format())) +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% unique(inputData_BioSamples[, Sample.Type]) ),
                            ggplot2::aes( x = -2, y = get(dependent),  shape = Sample.Type, colour = get(ClassCol)), size = 2, na.rm = TRUE) #+#, shape = 1, col = "purple"
      legend_order <- c(unique(inputData_BioSamples$Sample.Type))
    }
    #ifelse(!is.na(Class),get(Class), "black")

    if(!is.null(inputData_QC )){
      legend_order <- c(legend_order, unique(inputData_QC$Sample.Type))

      plotlinearData <-  plotlinearData +
        ggplot2::geom_point(data = subset(data_Signals, Sample.Type %in% QCs$Sample.Type),
                            ggplot2::aes( x = x, y = get(dependent),  shape = as.factor(Sample.Type),
                                          color = get(ClassCol)), size = 2, na.rm = TRUE) +
        ggplot2::scale_shape_manual(values = 1: (nlevels(as.factor(QCs$Sample.Type)) +1), breaks=rev(legend_order))

      #+, shape = 5

    }


    if(!is.null(TRANSFORM_Y)){
      plotlinearData <- plotlinearData +
        ggplot2:: scale_y_continuous(name = "Area", labels = function(l) {trans = scales::scientific(get(inverse_y)(l)); paste0(TRANSFORM_Y,"(", trans, ")")})#{paste0(TRANSFORM_Y,"(", scales::scientific(get(inverse_y)(l), ")"))})
      #labels = scales::trans_format(get(inverse_y)),
      #limits = c(NA, ggplot2::layer_scales(plotlinearData)$y$get_limits()[2] + 1 ))
      if("signalBlankRatio" %in% colnames(data_Signals)){

        plotlinearData <- plotlinearData +
          ggplot2::geom_hline(ggplot2::aes(yintercept = get(TRANSFORM_Y)(medBlank*signal_blank_ratio)), color = "grey", linewidth = 0.5)
        #   ggplot2::geom_hline(yintercept = get(TRANSFORM_Y)(unique(na.omit(data_Signals$medBlank))*signal_blank_ratio), color = "black") +
        #   ggplot2::annotate(x = max(data_Signals[[indipendent]], na.rm = T)/2 +2, y = get(TRANSFORM_Y)(data_Signals$medBlank*signal_blank_ratio)+0.2,geom = "text", label = paste(signal_blank_ratio," x median blank", size = 8, color = "black"))
      }

    } else{

      plotlinearData <- plotlinearData +
        ggplot2:: scale_y_continuous(name = "Area", labels = scales::scientific_format())

      if("signalBlankRatio" %in% colnames(data_Signals)){

        plotlinearData <- plotlinearData +
          ggplot2::geom_hline(ggplot2::aes(yintercept = medBlank*signal_blank_ratio), color = "grey", linewidth = 0.5)
        #   ggplot2::geom_hline(yintercept = get(TRANSFORM_Y)(unique(na.omit(data_Signals$medBlank))*signal_blank_ratio), color = "black") +
        #   ggplot2::annotate(x = max(data_Signals[[indipendent]], na.rm = T)/2 +2, y = get(TRANSFORM_Y)(data_Signals$medBlank*signal_blank_ratio)+0.2,geom = "text", label = paste(signal_blank_ratio," x median blank", size = 8, color = "black"))
      }

      # if("signalBlankRatio" %in% colnames(data_Signals)){
      #
      #   plotlinearData <- plotlinearData +
      #     ggplot2::geom_hline(yintercept = unique(na.omit(data_Signals$medBlank))*signal_blank_ratio, color = "black") +
      #     ggplot2::annotate(x = max(data_Signals[[indipendent]], na.rm = T)/2 +2, y = data_Signals$medBlank*signal_blank_ratio+200,geom = "text", label = paste(signal_blank_ratio," x median blank", size = 8, color = "black"))



    }




    if(length(GroupIndices) > 1 | any(GroupIndices %in% "all") | length(Feature) > 1 | any(Feature %in% "all" )){

      if(data.table::uniqueN(data_Signals[[Col_Batch]]) > 1){
        plotlinearData <-  plotlinearData +
          ggforce::facet_grid_paginate(ID ~ Batch ,
                                       scales = "free", ncol = nCol,nrow = nrRow, page = 1 )
      } else {
        plotlinearData <-  plotlinearData +
          ggforce::facet_grid_paginate(as.character(ID) ~. ,
                                       scales = "free", ncol = nCol,nrow = nrRow, page = 1 )
      }
    }


    if("signalBlankRatio" %in% colnames(data_Signals)){

      blankmedian <- data.frame(
        data_Signals |> dplyr::select(groupIndices, ID, Batch) |> unique()
      )


      blankmedian$label <- paste(signal_blank_ratio," x median blank")
      blankmedian$xtext <- data_Signals |>
        dplyr::group_by(groupIndices) |>
        dplyr::reframe(.groups = "drop", xblank = max(get(indipendent), na.rm = T)/2 + 2) |>
        dplyr::select(xblank) |>
        unlist(use.names = F)
      blankmedian <- data_Signals |>
        dplyr::group_by(groupIndices) |>
        dplyr::reframe(.groups = "drop",ytext = na.omit(medBlank*signal_blank_ratio)) |>
        unique() |>
        dplyr::ungroup() |>
        dplyr::right_join(blankmedian, by = "groupIndices")
      #dplyr::select(yblank) |>
      #unlist(use.names = F)


      if(!is.null(TRANSFORM_Y)){

        plotlinearData <- plotlinearData +

          ggplot2::geom_text(data = blankmedian, mapping = ggplot2::aes(x = xtext, y = get(TRANSFORM_Y)(ytext) + 0.5, label = label),size = 2) #+
      } else{
        ggplot2::geom_text(data = blankmedian, mapping = ggplot2::aes(x = xtext, y = ytext + 50, label = label),size = 2) #+
      }


    }



    if(printR2 %in% TRUE & any(!is.na(data_Signals$R2))) {

      text_label <- data_Signals[IsLinear %in% TRUE, .(ID, Batch, R2, y = get(dependent))]
      text_label <- text_label |> dplyr::group_by(ID, Batch) |>
        dplyr::mutate(R2 = ifelse(any(!is.na(R2)), R2[!is.na(R2)] , NA),
                      max_y = max(y, na.rm = T)) |>
        dplyr::select(-y) |>
        unique()




      plotlinearData <- plotlinearData +
        ggplot2::geom_text(data = text_label,
                           #ggplot2::geom_text(data = text_label[1:20,],#subset(data_Signals, !is.na(R2)),
                           ggplot2::aes(x = 0, y = max_y + 5, label = paste(Series,": R2 = ", round(R2,4))),
                           size = 3,
                           hjust = -0.1,
                           vjust = 2,
                           #inherit.aes = FALSE
        )


    }





    if(!is.null(inputData_BioSamples ) | !is.null(inputData_QC )){

      plotlinearData <-  plotlinearData +
        ggplot2::annotate(geom = "text", x = -3.5, y = Inf, label = "QC & Samples", size = 3,vjust = 2, na.rm = TRUE) #+
      #ggplot2::scale_color_manual(name = "In linear Range:",
      #                   values = c("FALSE" = "red", "TRUE" = "purple")) +
      #ggplot2::scale_shape_manual(values = c(0, 2, 5, 1, 6, 9, 3,4,7,8), breaks=rev(legend_order))
    }

    plotlinearData <-  plotlinearData +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
      ggplot2::theme(panel.grid.major=ggplot2::element_blank()) +
      ggplot2::theme(panel.background=ggplot2::element_blank()) +
      ggplot2::theme(axis.line=ggplot2::element_line()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90)) +
      ggplot2::theme(legend.position="top") +
      ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                      shape = ggplot2::guide_legend(order = 2))


    if(printPDF %in% TRUE){

      #plotObj <- vector("list", npage)
      if(data.table::uniqueN(data_Signal$Batch) <= 2){
        pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 9, height = 15)
        nrRow = 10
      }else {
        pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)
        nrRow = 5
      }

      n <- ggforce::n_pages(plotlinearData)

      for(i in 1:n) {

        if(data.table::uniqueN(data_Signal$Batch) > 1){
          suppressWarnings(print(plotlinearData + ggforce::facet_grid_paginate(ID ~ Batch ,
                                                                               scales = "free", ncol = nCol,nrow = nrRow, page = i )))
        } else {
          suppressWarnings( print(plotlinearData + ggforce::facet_grid_paginate(ID ~.,
                                                                                scales = "free", ncol = nCol,nrow = nrRow, page = i )))

        }
        print(paste0(i, " of ", n))



      }

      dev.off()

    }

    return(suppressWarnings(print(plotlinearData)))
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
  ID <- COLNAMES[["Feature_ID"]]
  Col_Batch = COLNAMES[["Batch"]]
  x <- X
  y <- Y

  data.table::setDT(inputData_Series)

  data_Signals_summary <- inputData_Series |>
    dplyr::group_by(DilutionPoint, Batch = get(Col_Batch)) |>
    dplyr::reframe(Missing = sum(is.na(get(y))),
              Blank = sum(signalBlankRatio, na.rm = TRUE),
              OutlierFOD = sum(OutlierFOD, na.rm = T),
              Trim = sum(trim %in% TRUE, na.rm = T),
              OutlierSOD = sum(OutlierSOD, na.rm = T),
              Slope = sum(trimPos %in% TRUE, na.rm = T),
              LR_TRUE = sum(IsLinear %in% TRUE, na.rm = T),
              LR_FALSE = sum(IsLinear %in% FALSE, na.rm = T)
              ) |>
    dplyr::group_by(DilutionPoint, Batch = get(Col_Batch)) |>
    dplyr::mutate(min_Feature = data.table::uniqueN(inputData_Series[[ID]]) - sum(Missing, Blank, OutlierFOD, Trim, OutlierSOD, Slope, LR_TRUE, LR_FALSE)) |>
    tidyr::pivot_longer(names_to = "Type",values_to =  "count" ,cols = -c(1:2))



  data_Signals_summary$Type <- factor(data_Signals_summary$Type, levels = rev(c("LR_TRUE", "LR_FALSE", "Slope", "OutlierSOD", "Trim", "OutlierFOD", "Blank", "Missing", "min_Feature")))
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
    ggplot2::scale_fill_manual(name = "", values = c('#822D3F','#8c510a','#d8b365','#A79EE5','#f6e8c3','#F194C1','#c7eae5','#5ab4ac','#01665e')) +
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
                                        COLNAMES, X, Y , output_dir,
                                        group = "Sample.Type",
                                        group2 = NULL,
                                        ordered = "Sample_ID"){


  assertthat::not_empty(inputData_Samples)
  # LR_object,statusLinear = c(TRUE, FALSE),
  ID <- COLNAMES[["Feature_ID"]]
  Col_Batch <- COLNAMES[["Batch"]]

  x <- X
  y <- Y
  SAMPLE_ID <- COLNAMES[["Sample_ID"]]
  ClassGroup <- COLNAMES[["Class"]]

  data.table::setDT(inputData_Samples)

  if(is.null(group2)){group2 <- group}

  if(group2 != group & sum(data.table::uniqueN(inputData_Samples[[group]]), data.table::uniqueN(inputData_Samples[[group2]])) > 4 ){
    nrRow <- 4
  } else if(group2 == group & data.table::uniqueN(inputData_Samples[[group]]) > 4 ){
    nrRow <- 4
  } else if (group2 == group) {
    nrRow = data.table::uniqueN(inputData_Samples[[group]])
  } else{
    nrRow = sum(data.table::uniqueN(inputData_Samples[[group]]), data.table::uniqueN(inputData_Samples[[group2]]))
  }

  inputData_Samples$ID = inputData_Samples[[ID]]
  inputData_Samples$Batch = inputData_Samplesl[[Col_Batch]]
  inputData_Samples$Sample.Type = inputData_Samples[[Sample.Type]]
  inputData_Samples$Batch = inputData_Samplesl[[Col_Batch]]
  inputData_Samples$Class = inputData_Samples[[ClassGroup]]
  inputData_Samples$Sample_ID = inputData_Samplesl[[SAMPLE_ID]]


  data_Signals_sample_summary <- inputData_Samples |>
    dplyr::group_by(Sample_ID, Batch, Sample.Type, Class, PlotOrder = get(ordered), group_2 = get(group2)) |>
    dplyr::reframe(
      Missing = sum(is.na(y), na.rm = T),
      LR_TRUE = sum(Status_LR %in% TRUE, na.rm = T),
      LR_FALSE = sum(Status_LR %in% c(FALSE, NA), na.rm = T),
      LR_ULOL = sum(Status_LR %in% "ULOL", na.rm = T),
      LR_BLOL = sum(Status_LR %in% "BLOL", na.rm = T),
      LR_TRUE_perc = LR_TRUE /(Missing + LR_TRUE + LR_FALSE + LR_ULOL + LR_BLOL),
      LR_ULOL_perc = LR_ULOL /(Missing + LR_TRUE + LR_FALSE + LR_TRUE + LR_BLOL),
      LR_BLOL_perc = LR_ULOL /(Missing + LR_TRUE + LR_FALSE + LR_TRUE + LR_ULOL)
    ) |>
    tidyr::pivot_longer(names_to = "Type",values_to =  "count" ,cols = -c(1:6))


  data_Signals_sample_summary$Type <- factor(data_Signals_sample_summary$Type, levels = rev(c("LR_TRUE_perc","LR_ULOL_perc","LR_BLOL_perc","LR_TRUE", "LR_FALSE","LR_ULOL" , "LR_BLOL", "Missing")))
  #data_Signals_sample_summary <- data_Signals_sample_summary |> dplyr::filter(count > 0)
  setorderv(data_Signals_sample_summary, "PlotOrder")
  data_Signals_sample_summary$Sample_ID <- factor(data_Signals_sample_summary$Sample_ID, levels = unique(data_Signals_sample_summary$Sample_ID))


  plot_Summary_samples <-
    ggplot2::ggplot(data = data_Signals_sample_summary|> dplyr::filter(Type %in% c("LR_TRUE_perc","LR_ULOL_perc" , "LR_BLOL_perc")),
                    ggplot2::aes(x = Sample_ID, y = count, label = count, fill = Type)) +
    ggplot2::geom_bar(stat="identity",
                      position="stack",
                      width = 0.5 )# +
  #geom_text(size = 3, position = position_fill(vjust = 0.5)) +
  if(group2 == group){ plot_Summary_samples <- plot_Summary_samples +
    ggforce::facet_grid_paginate(. ~ get(group), scales = "free", ncol = 1, nrow = nrRow, page = 1 )
  } else{
    plot_Summary_samples <- plot_Summary_samples +
      ggforce::facet_grid_paginate(group_2 ~ get(group), scales = "free", ncol = 1,nrow = nrRow, page = 1 )
  }



  plot_Summary_samples <- plot_Summary_samples +
    #scale_x_continuous(breaks = data_Signals$DilutionPoint) +
    ggplot2:: scale_y_continuous(labels = scales::percent, name = "No of Signals [%]", limits = c(0,NA)) +
    ggplot2::scale_fill_manual(name = "", values = c('#c7eae5','#5ab4ac','#01665e','#822D3F','#8c510a','#d8b365','#A79EE5','#f6e8c3','#F194C1')) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major.x =ggplot2::element_blank()) +
    ggplot2:: theme(panel.background=ggplot2::element_blank()) +
    ggplot2:: theme(axis.text.x = ggplot2::element_text(angle=90)) +
    ggplot2:: theme(axis.line=ggplot2::element_line()) #+
  #ggplot2::theme(legend.position = "none")

  if(printPDF %in% TRUE){

    pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)


    n <- ggforce::n_pages(plot_Summary_samples)

    for(i in 1:n) {

      if(group2 == group){
        suppressWarnings(print(plot_Summary_samples + ggforce::facet_grid_paginate(. ~ get(group) ,
                                                                             scales = "free", ncol = 1,nrow = nrRow, page = i )))
      } else {
        suppressWarnings( print(plot_Summary_samples + ggforce::facet_grid_paginate(group_2 ~ get(group),
                                                                              scales = "free", ncol = 1,nrow = nrRow, page = i )))

      }
      print(paste0(i, " of ", n))



    }

    dev.off()

  }



  return(  suppressWarnings(plot( plot_Summary_samples)))


}





#' Combine multiple omics signal tables into a single dataset
#'
#' This function merges different input data tables (series data, biological samples,
#' and QC samples) into a single unified data frame using common column intersections.
#' The merging is performed using full joins to retain all available observations.
#'
#' @param inputData_Series Data frame containing primary series-level measurements.
#' @param inputData_BioSamples Optional data frame containing biological sample data.
#' @param inputData_QC Optional data frame containing quality control (QC) sample data.
#'
#' @return A merged data frame containing all input datasets combined by shared columns.
#' @export
#'


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



#' Feature dilution series (FDS) diagnostic and calibration plotting
#'
#' This function generates detailed diagnostic plots for feature dilution series (FDS)
#' experiments. It visualizes dilution curves, linear range behavior, QC and biological
#' sample distributions, outlier annotations, and optional diagnostic plots (residuals
#' and QQ plots). The function can optionally export multi-page PDF reports.
#'
#' Key features:
#' - Visualization of dilution-dependent signal behavior per feature
#' - Highlighting of linear range, outliers, and QC performance
#' - Support for biological samples and QC overlays
#' - Optional diagnostic plots (residuals and QQ plots)
#' - Automatic faceting by feature and batch
#'
#' @param inputData_Series Data frame containing series dilution measurements and annotations.
#' @param inputData_BioSamples Optional data frame of biological sample measurements.
#' @param inputData_QC Optional data frame of QC sample measurements.
#' @param nrFeature Maximum number of features to include in the plot.
#' @param signal_blank_ratio Numeric; multiplier for blank signal threshold visualization.
#' @param printPDF Logical; if TRUE, plots are exported as a PDF file.
#' @param GroupIndices Character or vector; subset of group indices to plot ("all" for no filtering).
#' @param Feature Character or vector; subset of features to plot ("all" for no filtering).
#' @param printR2 Logical; if TRUE, R² and related statistics are displayed.
#' @param outputfileName Character; name of the output PDF file (without extension).
#' @param TRANSFORM_Y Function name used to transform y-values (e.g. "log", "exp").
#' @param inverse_y Function name used for inverse transformation of y-values.
#' @param COLNAMES Named list defining required column mappings (e.g. Feature_ID, Batch, Class, Sample_type).
#' @param Xcol Column name used as independent variable.
#' @param Ycol Column name used as dependent variable.
#' @param Series Optional label for the QC dilution series.
#' @param output_dir Directory where output PDF files are saved.
#' @param diagnostic Logical; if TRUE, residual and QQ diagnostic plots are generated.
#' @param ...
#' @return A patchwork ggplot object containing the dilution series visualization and optional diagnostics.
#' @export
#'

plot_FDS <- function(inputData_Series,
                     inputData_BioSamples,
                     inputData_QC,
                     nrFeature = 50,
                     signal_blank_ratio = 5,
                     printPDF = TRUE,
                     GroupIndices = "all",
                     Feature = "all",
                     printR2 = TRUE,
                     outputfileName = c("Calibrationplot"),
                     TRANSFORM_Y = "log",
                     inverse_y = "exp",
                     COLNAMES,
                     Xcol,
                     Ycol,
                     Series = "QC dilution Curves",
                     output_dir,
                     diagnostic = TRUE ){

  # ----------------------------
  # Input Validation
  # ----------------------------
  assertthat::not_empty(inputData_Series)

  # Spaltenzuordnung aus Mapping
  ID <- COLNAMES[["Feature_ID"]]
  Col_Batch <- COLNAMES[["Batch"]]
  independent <- Xcol
  dependent <- Ycol
  ClassCol <- COLNAMES[["Class"]]
  Sample.Type <- COLNAMES[["Sample_type"]]

  # ----------------------------
  # Converting in data.table
  # ----------------------------
  data.table::setDT(inputData_Series)
  inputData_Series$Sample.Type <- inputData_Series[[Sample.Type]]

  # ----------------------------
  # prepare BioSamples
  # ----------------------------
  if(!is.null(inputData_BioSamples)){
    inputData_BioSamples$Sample.Type <- inputData_BioSamples[[Sample.Type]]
    data.table::setDT(inputData_BioSamples)

    # Y-Transformation
    if(!is.null(TRANSFORM_Y)){
      inputData_BioSamples$Sample_area <-
        get(TRANSFORM_Y)(inputData_BioSamples[[COLNAMES[["Y"]]]])

      data.table::setnames(inputData_BioSamples,
                           old = "Sample_area",
                           new = dependent)

      inputData_BioSamples <- inputData_BioSamples |>
        dplyr::select(which(!duplicated(names(inputData_BioSamples))))
    }
  }

  # ----------------------------
  # prepare QC Data
  # ----------------------------
  if(!is.null(inputData_QC)){
    inputData_QC$Sample.Type <- inputData_QC[[Sample.Type]]
    data.table::setDT(inputData_QC)

    if(!is.null(TRANSFORM_Y)){
      inputData_QC$Sample_area <-
        get(TRANSFORM_Y)(inputData_QC[[COLNAMES[["Y"]]]])

      data.table::setnames(inputData_QC,
                           old = "Sample_area",
                           new = dependent)

      inputData_QC <- inputData_QC |>
        dplyr::select(which(!duplicated(names(inputData_QC))))
    }
  }

  # ----------------------------
  # combine Data
  # ----------------------------
  data_Signal <- combineData(inputData_Series,
                             inputData_BioSamples,
                             inputData_QC)

  # IDs setzen
  data_Signal$ID <- data_Signal[[ID]]
  data_Signal$Batch <- data_Signal[[Col_Batch]]
  data_Signal$Sample.Type <- data_Signal[[Sample.Type]]

  # ----------------------------
  # QC positioning on X-Axis
  # ----------------------------
  if(!is.null(inputData_QC)) {

    QCs <- data.frame(
      Sample.Type = unique(inputData_QC[[Sample.Type]]),
      x = -seq_along(unique(inputData_QC[[Sample.Type]]))
    )

    data_Signal <- dplyr::full_join(QCs, data_Signal, by = "Sample.Type")
  } else {
    data_Signal$x <- 0
  }

  # ----------------------------
  # sort
  # ----------------------------
  data.table::setorderv(data_Signal, c(ID, Col_Batch))

  # ----------------------------
  # Feature / Group Filter
  # ----------------------------
  if(any(GroupIndices != "all" & GroupIndices != "")){
    data_Signal <- data_Signal[groupIndices %in% GroupIndices]
  }

  if(any(Feature != "all" & Feature != "")){
    data_Signal <- data_Signal[get(ID) %in% Feature]
  }

  # ----------------------------
  # Feature Sampling
  # ----------------------------
  if(any(Feature %in% "all") &
     any(GroupIndices %in% "all") &
     data.table::uniqueN(data_Signal[[ID]]) > nrFeature){

    # ausgewogenes Sampling zwischen TRUE/FALSE LR Status
    ...

  }

  # ----------------------------
  # Plot Layout Parameter
  # ----------------------------
  nrRow <- ifelse(data.table::uniqueN(data_Signal$groupIndices) >= 4, 4,
                  data.table::uniqueN(data_Signal$groupIndices))

  nCol <- data.table::uniqueN(data_Signal[[Col_Batch]])

  # ----------------------------
  # color definitions
  # ----------------------------
  class_levels <- unique(stats::na.omit(data_Signal[[ClassCol]]))
  sample_levels <- unique(inputData_BioSamples[[Sample.Type]])
  batch_levels <- unique(data_Signal$Batch)

  class_colors <- stats::setNames(scales::hue_pal()(length(class_levels)), class_levels)
  sample_shapes <- stats::setNames(c(16,17,15,18,19)[seq_along(sample_levels)], sample_levels)
  batch_fills <- stats::setNames(RColorBrewer::brewer.pal(max(length(batch_levels),3),"Set1"),
                          batch_levels)

  # ----------------------------
  # main plot-Function pro Feature
  # ----------------------------
  metabolite_row <- function(df, show_x = FALSE, show_legend = FALSE, show_title = FALSE){

    # Einzel-Feature Plot (Dilution curve)
    plotlinearData <- ggplot2::ggplot(data = df,
                                      mapping = ggplot2::aes(x = get(independent), y = get(dependent))) +
      ggplot2::geom_point(na.rm = TRUE)

    # Outlier Markierungen
    if("OutlierFOD" %in% names(df)) {
      plotlinearData <- plotlinearData +
        ggplot2::geom_point(data = df[data = df$OutlierFOD == TRUE,], color="red", na.rm = TRUE)
    }

    # InRange + Fit Linie
    if("InRange" %in% names(df)){
      plotlinearData <- plotlinearData +
        ggplot2::geom_line(data = df[df$InRange == TRUE,], mapping = ggplot2::aes(y = predicted), na.rm = TRUE)
    }

    return(plotlinearData)
  }

  # ----------------------------
  # PDF Export
  # ----------------------------
  if(printPDF){

    grDevices::pdf(file.path(output_dir,
                  paste0(Sys.Date(),"_",outputfileName,".pdf")),
        width = 15, height = 9)

    metabolite_list <- split(data_Signal, data_Signal$ID)

    for(page in split(metabolite_list,
                      ceiling(seq_along(metabolite_list)/nrRow))){

      plots <- lapply(page, metabolite_row)

      print(patchwork::wrap_plots(plots))

    }

    grDevices::dev.off()
  }
}

#' Summary bar plot of feature-level quality and filtering status across dilution series
#'
#' This function generates a stacked bar plot summarizing different feature filtering
#' categories (e.g. missing values, blanks, outliers, linear range classification)
#' across dilution points and batches. Optionally, the plot can be exported as a PDF.
#'
#' @param inputData_Series Data frame containing series-level intensity and annotation data.
#' @param printPDF Logical; if TRUE, the plot is saved as a PDF file.
#' @param GroupIndices Character or vector; optional filtering of group indices ("all" for no filtering).
#' @param Feature Character or vector; optional selection of features ("all" for no filtering).
#' @param outputfileName Character; name of the output PDF file (without extension).
#' @param COLNAMES Named list defining column mappings (e.g. Feature_ID, Batch).
#' @param X Column name used as x-axis variable (e.g. dilution index).
#' @param Y Column name used as y-axis variable (e.g. signal intensity).
#' @param output_dir Directory where the PDF file will be saved.
#'
#' @return A ggplot object representing the summary bar plot.
#' @export

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
              Blank = ifelse("signalBlankRatio" %in% colnames(inputData_Series), sum(signalBlankRatio, na.rm = TRUE), 0),
              OutlierFOD =  ifelse("OutlierFOD" %in% colnames(inputData_Series), sum(OutlierFOD, na.rm = T), 0),
              Trim = ifelse("trim" %in% colnames(inputData_Series), sum(trim %in% TRUE, na.rm = T), 0),
              OutlierSOD = ifelse("OutlierSOD" %in% colnames(inputData_Series), sum(OutlierSOD, na.rm = T), 0),
              #Slope = sum(trimPos %in% TRUE, na.rm = T),
              LR_TRUE = sum(InRange %in% TRUE, na.rm = T),
              LR_FALSE = sum(InRange %in% FALSE, na.rm = T)
              ) |>
    dplyr::group_by(DilutionPoint, Batch = get(Col_Batch)) |>
    dplyr::mutate(min_Feature = data.table::uniqueN(inputData_Series[[ID]]) - sum(Missing, Blank, OutlierFOD, Trim, OutlierSOD,  LR_TRUE, LR_FALSE)) |> #Slope,
    tidyr::pivot_longer(names_to = "Type",values_to =  "count" ,cols = -c(1:2))



  data_Signals_summary$Type <- factor(data_Signals_summary$Type, levels = rev(c("LR_TRUE", "LR_FALSE", "OutlierSOD", "Trim", "OutlierFOD", "Blank", "Missing", "min_Feature"))) #"Slope",
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

    grDevices::pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)}


  print( plot_Summary)

  if(printPDF %in% TRUE){ grDevices::dev.off()}

  return( plot_Summary)

}

#' Summary bar plot of signal quality per sample
#'
#' This function generates a stacked bar plot showing the proportion and counts
#' of different signal quality categories (e.g. linear range status, missing values,
#' out-of-linear-range signals) per sample. The results can be grouped by sample type
#' and optionally a second grouping variable, and exported as a multi-page PDF.
#'
#' @param inputData_Samples Data frame containing sample-level annotation and signal data.
#' @param printPDF Logical; if TRUE, the plot is saved as a PDF file.
#' @param GroupIndices Character or vector; optional filtering of group indices ("all" for no filtering).
#' @param Feature Character or vector; optional selection of features ("all" for no filtering).
#' @param outputfileName Character; name of the output PDF file (without extension).
#' @param COLNAMES Named list defining column mappings (e.g. Feature_ID, Batch, Sample_type, Sample_ID, Class).
#' @param X Column name used as x-axis variable.
#' @param Y Column name used as y-axis variable.
#' @param output_dir Directory where the PDF file will be saved.
#' @param group Primary grouping variable used for faceting (default: "Sample.Type").
#' @param group2 Optional secondary grouping variable; if NULL, defaults to `group`.
#' @param ordered Column name used to define plotting order of samples.
#'
#' @return A ggplot object of the sample summary bar plot.
#' @export



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
  Sample.Type <- COLNAMES[["Sample_type"]]
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
  inputData_Samples$Batch = inputData_Samples[[Col_Batch]]
  inputData_Samples$Sample.Type = inputData_Samples[[Sample.Type]]
  inputData_Samples$Batch = inputData_Samples[[Col_Batch]]
  inputData_Samples$Class = inputData_Samples[[ClassGroup]]
  inputData_Samples$Sample_ID = inputData_Samples[[SAMPLE_ID]]


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
  data.table::setorderv(data_Signals_sample_summary, "PlotOrder")
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

    grDevices::pdf(file = file.path(output_dir,paste0(Sys.Date(),"_", outputfileName,".pdf")), width = 15, height = 9)


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

    grDevices::dev.off()

  }



  return(suppressWarnings(plot( plot_Summary_samples)))


}

#' Title
#'
#' @param inputData
#' @param Compound_ID
#' @param Sample_ID
#' @param Y
#' @param inputDataSummary
#' @param outputformat
#' @param PATH
#' @param filter.type
#' @param filter.calc
#' @param filter.calc.x.Batch
#' @param filter.calc.x.Class
#' @param filter.percent
#'
#' @return
#' @export
#'
#' @examples
MS_filterSamples <- function(
    inputData,
    Compound_ID,
    Sample_ID,
    Y,
    inputDataSummary,
    outputformat = c("wide", "long", "both"),
    PATH,
    filter.type = c("Batch", "Class", "Batch_Class"),
    filter.calc = c("all", "at least one", "at least x"),
    filter.calc.x.Batch = NULL,
    filter.calc.x.Class = NULL,
    filter.percent = 80
){

  data.table::setDT(inputData)
  data.table::setDT(inputDataSummary)

  dat <- data.table::copy(inputData)
  datSummary <- data.table::copy(inputDataSummary)

  datSummary <- datSummary[!grep("notEnoug", LRFlag)]


  if(filter.type %in% "Batch"){

    BatchSummary <- datSummary |>
      dplyr::group_by(ID, Batch) |>
      dplyr::summarize(.groups = "drop",
                       LR_TRUE_Batch = round(sum(LR_TRUE, na.rm = TRUE)/sum(Signals)*100,2))

    if(filter.calc %in% "all") {

      Compounds <- BatchSummary |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter = all(LR_TRUE_Batch >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(ID) |> unique()
    }else if(filter.calc %in% "at least one"){
      Compounds <- BatchSummary |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter = any(LR_TRUE_Batch >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(ID) |> unique()
    }else if(filter.calc %in% "at least x"){
      Compounds <- BatchSummary |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter = sum(LR_TRUE_Batch >= filter.percent)) |>
        dplyr::filter(Filter >= filter.calc.x.Batch) |>
        dplyr::select(ID) |> unique()
    }
  }else if(filter.type %in% "Class"){

    ClassSummary <- datSummary |>
      dplyr::group_by(ID, Class) |>
      dplyr::summarize(.groups = "drop",
                       LR_TRUE_Class = round(sum(LR_TRUE, na.rm = TRUE)/sum(Signals)*100,2))

    if(filter.calc %in% "all") {

      Compounds <- ClassSummary |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter = all(LR_TRUE_Class >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(ID) |> unique()
    }else if(filter.calc %in% "at least one"){
      Compounds <- ClassSummary |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter = any(LR_TRUE_Class >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(ID) |> unique()
    }else if(filter.calc %in% "at least x"){
      Compounds <- ClassSummary |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter = sum(LR_TRUE_Class >= filter.percent)) |>
        dplyr::filter(Filter >= filter.calc.x.Batch) |>
        dplyr::select(ID) |> unique()
    }
  }else if(filter.type %in% "Batch_Class"){

    BatchClassSummary <- datSummary |>
      dplyr::group_by(ID, Batch,Class) |>
      dplyr::summarize(.groups = "drop",
                       LR_TRUE_BatchClass = round(sum(LR_TRUE, na.rm = TRUE)/sum(Signals)*100,2))

    if(filter.calc %in% "all") {
      # all classes in at least x Batches
      Compounds <- BatchClassSummary |>
        dplyr::group_by(ID, Batch) |>
        dplyr::mutate(Filter = all(LR_TRUE_BatchClass >= filter.percent)) |>
        dplyr::ungroup() |>
        dplyr::select(ID, Batch, Filter) |> unique() |>
        dplyr::group_by(ID) |>
        dplyr::mutate(Filter2 = sum(Filter, na.rm = TRUE) >= filter.calc.x.Batch) |>
        dplyr::filter(Filter2 %in% TRUE) |>
        dplyr::select(ID) |> unique()
    }else if(filter.calc %in% "at least one"){
      # at least one identical class in at least x Batches
      Compounds <- BatchClassSummary |>
        dplyr::group_by(ID, Class) |>
        dplyr::mutate(Filter = sum(LR_TRUE_BatchClass >= filter.percent)) |>
        dplyr::ungroup() |>
        dplyr::select(ID, Batch, Filter) |> unique() |>
        dplyr::group_by(ID, Batch) |>
        dplyr::mutate(Filter2 = any(Filter >= filter.calc.x.Batch)) |>
        dplyr::filter(Filter2 %in% TRUE) |>
        dplyr::ungroup() |>
        dplyr::select(ID) |> unique()
    }else if(filter.calc %in% "at least x"){
      # at least x identical classes in at least x Batches
      Compounds <- BatchClassSummary |>
        dplyr::group_by(ID, Class) |>
        dplyr::mutate(Filter = sum(LR_TRUE_BatchClass >= filter.percent) >= filter.calc.x.Batch) |>
        dplyr::ungroup() |>
        dplyr::group_by(ID, Batch) |>
        dplyr::mutate(Filter2 = sum(Filter) >= filter.calc.x.Class) |>
        dplyr::filter(Filter2 %in% TRUE) |>
        dplyr::ungroup() |>
        dplyr::select(ID) |> unique()
    }
  }

  if(outputformat %in% "long"){

    filteredData <- dat |> dplyr::filter(get(Compound_ID) %in% unlist(Compounds))

    writexl::write_xlsx(x = list(input = dat, filtered.long = filteredData),
                        path = file.path(PATH, paste0(Sys.Date(),"_filteredData.xlsx")))

    } else if(outputformat %in% "wide"){

      filteredData <- dat |> dplyr::filter(get(Compound_ID) %in% unlist(Compounds)) |>
        dplyr::select(ID = all_of(Compound_ID), Sample = all_of(Sample_ID), y = all_of(Y)) |>
        tidyr::pivot_wider(names_from = Sample, values_from = y)

      writexl::write_xlsx(x = list(input = dat, filtered.wide = filteredData),
                          path = file.path(PATH, paste0(Sys.Date(),"_filteredData.xlsx")))

  } else if(outputformat %in% "both"){

    filteredData <- list(
      wide = dat |> dplyr::filter(get(Compound_ID) %in% unlist(Compounds)) |>
      dplyr::select(ID = all_of(Compound_ID), Sample = all_of(Sample_ID), y = all_of(Y)) |>
      tidyr::pivot_wider(names_from = Sample, values_from = y),
       long = dat |> dplyr::filter(get(Compound_ID) %in% unlist(Compounds)))

    writexl::write_xlsx(x = list(input = dat,
                                 filtered.long = filteredData$long,
                                 filtered.wide = filteredData$wide),
                        path = file.path(PATH, paste0(Sys.Date(),"_filteredData.xlsx")))

  }


  return(filteredData)

}












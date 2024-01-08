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
    Compound_ID = "groupIndices",
    Sample_ID = "SampleID",
    col_Batch = "Batch",
    col_Class = "Group",
    Y,
    #inputDataSummary,
    outputformat = c("none", "wide", "long", "both"),
    PATH,
    filter.type = c("Batch", "Class", "Batch_Class"),
    filter.calc = c("all", "at least one", "at least x", "linear or 2 levels different"),
    filter.calc.x.Batch = NULL,
    filter.calc.x.Class = NULL,
    filter.percent = 80
){

  data.table::setDT(inputData)
  #data.table::setDT(inputDataSummary)

  dat <- data.table::copy(inputData)
  dat <- dat[!grep("notEnough", LRFlag)]


  #datSummary <- data.table::copy(inputDataSummary)

  #datSummary <- datSummary[!grep("notEnough", LRFlag)]


  if(filter.type %in% "Batch"){

    BatchSummary <- dat |>
      dplyr::group_by(Compound_ID = get(Compound_ID), col_Batch = get(col_Batch)) |>
      dplyr::reframe(  'Signals' = length(Status_LR),
                       'LR_Status' = vctrs::vec_count(Status_LR)$key,
                       'LR_Status_n' = vctrs::vec_count(Status_LR)$count,
                       'LR_Status_n[%]' = round(vctrs::vec_count(Status_LR)$count/length(Status_LR)*100,2))
    #LR_TRUE_Batch = round(sum(LR_TRUE, na.rm = TRUE)/sum(Signals)*100,2))

    if(filter.calc %in% "all") {

      Compounds <- BatchSummary |>
        dplyr::group_by(Compound_ID, LR_Status) |>
        dplyr::filter(LR_Status %in% c("TRUE", "ULOL", "BLOL")) |>
        dplyr::mutate(Filter =all(`LR_Status_n[%]` >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(Compound_ID) |> unique()
    }else if(filter.calc %in% "at least one"){
      Compounds <- BatchSummary |>
        dplyr::group_by(Compound_ID, LR_Status) |>
        dplyr::mutate(Filter = any(`LR_Status_n[%]` >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(Compound_ID) |> unique()
    }else if(filter.calc %in% "at least x"){
      Compounds <- BatchSummary |>
        dplyr::group_by(Compound_ID, LR_Status) |>
        dplyr::mutate(Filter = sum(`LR_Status_n[%]` >= filter.percent)) |>
        dplyr::filter(Filter >= filter.calc.x.Batch) |>
        dplyr::select(Compound_ID) |> unique()
    }
  }else if(filter.type %in% "Class"){

    ClassSummary <- dat |>
      dplyr::group_by(Compound_ID = get(Compound_ID), col_Class = get(col_Class)) |>
      dplyr::reframe(  'Signals' = length(Status_LR),
                       'LR_Status' = vctrs::vec_count(Status_LR)$key,
                       'LR_Status_n' = vctrs::vec_count(Status_LR)$count,
                       'LR_Status_n[%]' = round(vctrs::vec_count(Status_LR)$count/length(Status_LR)*100,2))


    if(filter.calc %in% "all") {

      Compounds <- ClassSummary |>
        dplyr::group_by(Compound_ID) |>
        dplyr::filter(LR_Status %in% "TRUE") |>
        dplyr::mutate(Filter = all(`LR_Status_n[%]` >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(Compound_ID) |> unique()

    }else if(filter.calc %in% "at least one"){

      Compounds <- ClassSummary |>
        dplyr::group_by(Compound_ID) |>
        dplyr::filter(LR_Status %in% "TRUE") |>
        dplyr::mutate(Filter = any(`LR_Status_n[%]` >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(Compound_ID) |> unique()



    }else if(filter.calc %in% "linear or 2 levels different"){
      Compounds_TRUE <- ClassSummary |>
        dplyr::group_by(Compound_ID) |>
        dplyr::filter(LR_Status %in% "TRUE") |>
        dplyr::mutate(Filter = any(`LR_Status_n[%]` >= filter.percent)) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(Compound_ID) |> unique()

      Compounds_outside <- ClassSummary |>
        dplyr::group_by(Compound_ID, LR_Status) |>
        dplyr::filter(!Compound_ID %in% Compounds_TRUE$Compound_ID) |>
        dplyr::filter(LR_Status %in% c("ULOL", "BLOL"),`LR_Status_n[%]` >= filter.percent) |>
        dplyr::group_by(Compound_ID, col_Class) |>
        dplyr::filter(`LR_Status_n[%]` == max(`LR_Status_n[%]`)) |>
        dplyr::group_by(Compound_ID) |>
        dplyr::mutate(Filter = data.table::uniqueN(col_Class) >= 2 & data.table::uniqueN(LR_Status) >= 2) |>
        dplyr::filter(Filter %in% TRUE) |>
        dplyr::select(Compound_ID) |> unique()

      Compounds <- rbind(Compounds_TRUE, Compounds_outside)

    }

    # }else if(filter.type %in% "Batch_Class"){
    #
    #   BatchClassSummary <- datSummary |>
    #     dplyr::group_by(ID, Batch,Class) |>
    #     dplyr::summarize(.groups = "drop",
    #                      LR_TRUE_BatchClass = round(sum(LR_TRUE, na.rm = TRUE)/sum(Signals)*100,2))
    #
    #   if(filter.calc %in% "all") {
    #     # all classes in at least x Batches
    #     Compounds <- BatchClassSummary |>
    #       dplyr::group_by(ID, Batch) |>
    #       dplyr::mutate(Filter = all(LR_TRUE_BatchClass >= filter.percent)) |>
    #       dplyr::ungroup() |>
    #       dplyr::select(ID, Batch, Filter) |> unique() |>
    #       dplyr::group_by(ID) |>
    #       dplyr::mutate(Filter2 = sum(Filter, na.rm = TRUE) >= filter.calc.x.Batch) |>
    #       dplyr::filter(Filter2 %in% TRUE) |>
    #       dplyr::select(ID) |> unique()
    #   }else if(filter.calc %in% "at least one"){
    #     # at least one identical class in at least x Batches
    #     Compounds <- BatchClassSummary |>
    #       dplyr::group_by(ID, Class) |>
    #       dplyr::mutate(Filter = sum(LR_TRUE_BatchClass >= filter.percent)) |>
    #       dplyr::ungroup() |>
    #       dplyr::select(ID, Batch, Filter) |> unique() |>
    #       dplyr::group_by(ID, Batch) |>
    #       dplyr::mutate(Filter2 = any(Filter >= filter.calc.x.Batch)) |>
    #       dplyr::filter(Filter2 %in% TRUE) |>
    #       dplyr::ungroup() |>
    #       dplyr::select(ID) |> unique()
    #   }else if(filter.calc %in% "at least x"){
    #     # at least x identical classes in at least x Batches
    #     Compounds <- BatchClassSummary |>
    #       dplyr::group_by(ID, Class) |>
    #       dplyr::mutate(Filter = sum(LR_TRUE_BatchClass >= filter.percent) >= filter.calc.x.Batch) |>
    #       dplyr::ungroup() |>
    #       dplyr::group_by(ID, Batch) |>
    #       dplyr::mutate(Filter2 = sum(Filter) >= filter.calc.x.Class) |>
    #       dplyr::filter(Filter2 %in% TRUE) |>
    #       dplyr::ungroup() |>
    #       dplyr::select(ID) |> unique()
    #   }
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

  filteredData <- dat |> dplyr::filter(get(Compound_ID) %in% unlist(Compounds))

  return(filteredData)

}



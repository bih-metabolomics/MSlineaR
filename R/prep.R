#' Prepare raw data
#'
#' @param dat Data frame with raw data. Mandatory columns are:
#' - ID,
#' - IDintern
#' - Batch,
#' - Replicate,
#' - MZ
#' - RT
#' - Concentration
#' - Intensity
#'
#'
#' @return Data frame with columns:
#' - IDintern,
#' - ID,
#' - IntensityNorm,
#' - DP,
#' - Comment,
#' - pch,
#' - color,
#' - ConcentrationRaw,
#' - IntensityRaw,
#' - ConcentrationLog,
#' - IntensityLog,
#' - IntensityRaw
#' - DilutionPoint
#'
#' @export
#'
#' @examples
#'
prepareData <- function(dat){

    dat <- dat %>%
    dplyr::rename(ConcentrationRaw = Concentration, IntensityRaw = Intensity) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(ConcentrationLog = log(ConcentrationRaw),
           DilutionPoint = row_number(),
           Comment = NA,
           IntensityNorm = IntensityRaw/max(IntensityRaw, na.rm = T)*100,
           IntensityLog = log(IntensityRaw),
           DP = DilutionPoint - mean(DilutionPoint)) %>%
    arrange(DilutionPoint) %>%
    ungroup() %>%
    mutate(pch = 19, color = "black") %>%
    dplyr::select(IDintern,
           ID,
           IntensityNorm,
           DP,
           Comment,
           pch,
           color,
           ConcentrationRaw,
           IntensityRaw,
           ConcentrationLog,
           IntensityLog,
           DilutionPoint)

  return(dat)
}

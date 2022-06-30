prepareData <- function(dat){

  dat <- dat %>%
    rename(ConcentrationRaw = Concentration, IntensityRaw = Intensity) %>%
    dplyr::group_by(ID) %>%
    mutate(ConcentrationLog = log(ConcentrationRaw),
           DilutionPoint = row_number(),
           Comment = NA,
           IntensityNorm = IntensityRaw/max(IntensityRaw, na.rm = T)*100,
           DP = DilutionPoint - mean(DilutionPoint)) %>%
    arrange(DilutionPoint) %>%
    ungroup() %>%
    mutate(pch = 19, color = "black") %>%
    select(IDintern,
           ID,
           IntensityNorm,
           DP,
           Comment,
           pch,
           color,
           ConcentrationRaw,
           IntensityRaw,
           ConcentrationLog,
           IntensityRaw,
           DilutionPoint)

  return(dat)
}

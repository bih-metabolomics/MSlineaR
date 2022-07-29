#' Title
#'
#' @param datList
#'
#' @return
#' @export
#'
#' @examples
func_hist_QC_summary <- function(datList = processList$dataLinearRange){

  NAdat <- datList$DilutionPoint[is.na(datList$IntensityRaw)] %>% table()
  NANr  <- lapply(1:length(unique(datList$DilutionPoint)), function(x) {
    ifelse(x %in% names(NAdat),yes = NAdat[names(NAdat) == x],no = 0)
  }) %>%  unlist()

  Outdat <- datList$DilutionPoint[datList$Outlier == TRUE] %>% table()
  OutNr  <- lapply(1:length(unique(datList$DilutionPoint)), function(x) {
    ifelse(x %in% names(Outdat),yes = Outdat[names(Outdat) == x],no = 0)
  }) %>%  unlist()

  ValidDat <- datList$DilutionPoint[datList$linearRange == TRUE] %>% table()
  ValidNr  <- lapply(1:length(unique(datList$DilutionPoint)), function(x) {
    ifelse(x %in% names(ValidDat),yes = ValidDat[names(ValidDat) == x],no = 0)
  }) %>%  unlist()


  library(gridExtra)

  summary_dat <- tibble(
    "Dilution Point" = 1 : length(unique(datList$DilutionPoint)),
    "Concentration" = unique(datList$ConcentrationRaw),
    "Nr of points within linear range" = ValidNr,
    "Nr of points without linear range" = rep(length(unique(datList$ID)), length(unique(datList$DilutionPoint))) - ValidNr,
    "Nr of NAs" = NANr,
    "Nr of Outliers" = OutNr,
  )

  g <-  summary_dat %>%
    gather(key = "Nr", value = "Values", 3:6) %>%
    #mutate(DP = `Dilution Points`) %>%

    ggplot( aes(x = `Dilution Point` , y = Values, fill = Nr )) +
    geom_bar(width = 0.5, position = position_dodge(), stat = "identity") +
    geom_hline(yintercept = n_distinct(datList %>% filter(`DilutionPoint` == 1)), linetype = "dotted") +
    theme_bw() +
    scale_x_continuous(breaks = seq(1, n_distinct(summary_dat$`Dilution Point`), by = 1),labels = unique(summary_dat$Concentration), name = "Dilution") #+
  #scale_y_continuous(breaks = seq(0, 12, by = 1))
  # annotate(geom = "table",
  #        x = 10,
  #        y = 0,
  #        label = list(summary_dat))
  tbl <- tableGrob(summary_dat, rows=NULL)

  grid.arrange(g, tbl,
               ncol = 1,
               as.table = TRUE)

}

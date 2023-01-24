
# MSlineaR

<!-- badges: start -->
<!-- badges: end -->

The goal of MSlineaR is to ...

## Installation

You can install the development version of MSlineaR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bih-metabolomics/MSlineaR")
```

## Example

Finding the linear Range for an targeted tryptophan Calibration data set:

``` r
library(MSlineaR)
library(data.table)

MSData <- MSlineaR::TryptophanData
Calibrants <- MSlineaR::TryptophanCalibration


data_tbl_cal <- data.table::data.table(MSData)[Sample.Type %in% "Calibration Standard", ]

data_tbl_cal <- data_tbl_cal[Calibrants, on = c(Sample.Identification = "Sample Identification", Compound = "Metabolite" )]
data_tbl_cal <- na.omit(data_tbl_cal, c("ConcentrationCal", "Batch"))

targetedMstarsCal <- AssessLinearity(

  DAT = data_tbl_cal,
  COLNAMES = c(
    ID =  "Compound",
    REPLICATE = "Batch",
    X = "ConcentrationCal",
    Y =   "Area"
  ),
  nCORE = 4,
  MIN_FEATURE = 5,
  LOG_TRANSFORM = TRUE,
  R2min = 0.9,
  res = 10
)




```


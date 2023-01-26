
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

# read in example data
MSData <- MSlineaR::TryptophanData
Calibrants <- MSlineaR::TryptophanCalibration

# select samples used for concentration or dilution series
data_tbl_cal <- data.table::data.table(MSData)[Sample.Type %in% "Calibration Standard", ]

# combine samples and concentrations
data_tbl_cal <- data_tbl_cal[Calibrants, on = c(Sample.Identification = "Sample Identification", Compound = "Metabolite" )]
data_tbl_cal <- na.omit(data_tbl_cal, c("Concentration", "Batch"))


targetedMSCal <- AssessLinearity(

  DAT = data_tbl_cal,
  COLNAMES = c(
    ID =  "Compound",
    REPLICATE = "Batch",
    X = "Concentration",
    Y =   "Area"
  ),
  nCORE = 4,
  MIN_FEATURE = 5,
  LOG_TRANSFORM = TRUE,
  R2min = 0.9,
  res = 10
)

# select samples
data_tbl_Sample = MSData[Sample.Type %in% "Sample"]

# check if samples within linear range or not
data_tbl_Sample = getLRstatus(dats = data_tbl_Sample, datCal = targetedMSCal$summaryFDS, COLNAMES = c(ID = "Compound", Replicate = "Batch", Y = "Area"))

# plot all FDS
plotFDS(LR_object = targetedMSCal,
        groupIndices = "all", 
        statusLinear = c(TRUE, FALSE), 
        ID = "all", 
        printPDF = TRUE, 
        outputfileName = c("Calibrationplot"), 
        pdfwidth = 9, 
        pdfheight = 60)
        
        
# plot all FDS plus Samples

plotSamples(LR_object = targetedMSCal,
            data_Sample = data_tbl_Sample,
            groupIndices = "all", 
            statusLinear = c(TRUE, FALSE), 
            ID = "all", 
            printPDF = TRUE, 
            outputfileName = c("Sampleplot"), 
            pdfwidth = 9, 
            pdfheight = 60)



```


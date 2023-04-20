# MSlineaR

<!-- badges: start -->

<!-- badges: end -->

The goal of `MSlineaR` is to clean up a data set of serial concentrated or diluted Compounds to only signals which show a linear response. It can be used for targeted and untargeted metabolomic data sets. In short the package use a serial diluted or concentrated data set to determine the linear portion of the dilution/calibration curve. Therefore the following steps will be performed:

1.  The data will be checked and transformed.

2.  Outliers will be detected using three curve fitting models.

3.  The beginning and end of dilution/calibration curve will be trimmed, so that the first point of the potential linear curve has the lowest intensity and the last point potential linear curve should have the largest intensity.

4.  A second outlier detection will be performed using the three curve fitting models.

5.  Reducing the curve to only consecutive points which showing a positive slope.

6.  Determining the linear part of the curve by using a linear regression model. Only curves with an Coefficient of Determination (R²) higher or equal to a predefined R² are considered as linear.

7.  If the data set consists of several batches, the algorithm checks if a linear part was found in all batches. If not than the compound will be flagged and if requested can be removed.

After these steps for each dilution/concentration curve are the information available about if the respective compound has a linear part and and where it starts and ends. With these information the signals of biological samples can be classified into signals which are in the linear range and signals which are outside of the linear range.

## Installation

To get the latest development version from `MSlineaR` from [GitHub](https://github.com/) you can use the [`devtools`](https://github.com/r-lib/devtools) package:

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

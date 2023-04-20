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

`MSlineaR` requires GNU R version \>= 4.1.2 installed.

To get the latest development version from `MSlineaR` from [GitHub](https://github.com/) you can use the [`devtools`](https://github.com/r-lib/devtools) package:

``` r
# install.packages("devtools")
devtools::install_github("bih-metabolomics/MSlineaR")
```

As long as `MSlineaR` is not a public package, you need to have access to my repository and additionally you need create a personal access token in github. Ask me if you have trouble with installing the package.

Install `MSlineaR` by running the following commands in R:

    if (!requireNamespace("usethis")) install.packages("usethis")
      
    # set config
    # add github user name and email adress, which you are using for github
    usethis::use_git_config(user.name = "USERNAME", user.email = "EMAIL")

    #Go to github page to generate token only for first time
    usethis::create_github_token()

    #paste your PAT into pop-up that follows...
    credentials::set_github_pat()

    # install.packages("devtools")
    devtools::install_github("bih-metabolomics/MSlineaR")

## Input Data

The input data needs at least 5 columns, indicating:

-   the Identification of the samples

-   the type of the samples, e.g. Calibrants, samples, pooled QC etc

-   the Identification of the compounds

-   the Instrument response( y) , e.g. Area, Intensity etc

-   the dilution step (1, 2, 3...) of a pooled QC for untargeted analysis and the Concentration of Calibrants for targeted analysis

| Sample_ID |  Type   | Compound_ID |  Area  | Dilution |
|:---------:|:-------:|:-----------:|:------:|:--------:|
|   cal1    | QC pool |  compound1  | 13578  |    1     |
|   cal2    | QC pool |  compound1  | 157898 |    2     |
|    s25    | sample  |  compound1  | 136789 |    NA    |

Optional the data table could include:

-   a column to distinguish between different Batches

-   a column to distinguish between different statistical classes e.g. healthy and treatment

## Example

Finding the linear Range for an targeted tryptophan Calibration data set:

``` r
library(MSlineaR)

data_tbl <- MSlineaR::targeted_MS


targetedMSCal <- AssessLinearity(
  dilution_factor = 3,
  analysis_type ="targeted",
  input_data = data_tbl,
  column_sample_type = "Sample.Type" ,
  sample_type_QC = c("Pooled QC","Reference QC","Blank"),
  sample_type_sample = "Sample",
  sample_type_serial = "Calibration Standard",
  sample_ID = "Sample.Identification",

  column_ID = "Compound",
  column_Batch = "Batch",
  column_X = "Concentration",
  column_Y = "Area",
  column_Y_sample = "Area",
  min_feature = 5,
  output_name = "Test_targeted",
  output_dir = "S:/Projects/2022/MSTARS/Janine_linearity_project/Benchmark data/targeted/230419",
  nCORE = 4

)
```

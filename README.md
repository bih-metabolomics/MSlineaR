# MSlineaR

<!-- badges: start -->

<!-- badges: end -->

The goal of `MSlineaR` is to clean up a data set of serial concentrated or diluted Compounds to only signals which show a linear response. It can be used for targeted and untargeted metabolomic data sets. In short the package use a serial diluted or concentrated data set to determine the linear portion of the dilution/calibration curve. Therefore the following steps will be performed:

1.  The data will be checked and transformed.

2.  The minimal lower limit of quantification will be determined using a blank and a user defined signal to blank ratio.

3.  Outliers will be detected using three curve fitting models.

4.  The beginning and end of dilution/calibration curve will be trimmed, so that the first point of the potential linear curve has the lowest intensity and the last point potential linear curve should have the largest intensity.

5.  A second outlier detection will be performed using the three curve fitting models.

6.  Reducing the curve to only consecutive points which showing a positive slope.

7.  Determining the linear part of the curve by using a linear regression model. Only curves with an Coefficient of Determination (R²) higher or equal to a predefined R² are considered as linear.

8.  If the data set consists of several batches, the algorithm checks if a linear part was found in all batches. If not than the compound will be flagged and if requested can be removed.

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

## Input Data Table

The input data needs at least 5 columns, indicating:

-   the Identification of the samples

-   the type of the samples, e.g. Calibrants, samples, pooled QC etc

-   the Identification of the compounds

-   the Instrument response( y) , e.g. Area, Intensity etc

-   the dilution step (1, 2, 3...) of a pooled QC for untargeted analysis and of the Concentration of Calibrants for targeted analysis

| Sample_ID |  Type   | Compound_ID |  Area  | Dilution |
|:---------:|:-------:|:-----------:|:------:|:--------:|
|   cal1    | QC pool |  compound1  | 13578  |    1     |
|   cal2    | QC pool |  compound1  | 157898 |    2     |
|    s25    | sample  |  compound1  | 136789 |    NA    |

Optional the input data table could include:

-   a column to distinguish between different Batches if there are more than one batch

-   a column to distinguish between different statistical classes e.g. healthy and treatment

Additional columns will be ignored for calculations, but will be present in the final output.

## Mandatory Input Arguments

-   *analysis_type:* You can choose between "targeted" or "untargeted"

-   *input_data*: data table or data frame with at least five columns, see above

    -   *column_sample_type*: Column in *input_data* which indicates the sample types of the data set, as:

        -   *sample_type_serial:* type name in *column_sample_type* of the serial diluted or concentrated samples

        -   *sample_type_blank*: type name in *column_sample_type* of the blank sample, which will be used to determine the minimum lower limit of quantification

    -   *column_sample_ID:* column in *input_data* indicating the sample ID

    -   *column_X:* column in *input_data* indicating the dilution or concentration of the serial diluted samples

    -   *column_Y:* column in *input_dat*a indicating the Instrument response, e.g. Area, Intensity etc

-   *dilution_factor*: distance between the dilution/concentration steps, e.g. dilution_factor = 2 means, that the concentration of dilution 1 is half the concentration of dilution 2

-   *signal_blank_ratio*: Default set to 5, according to EMA guidelines 2012, the signal should be at least 5 times the the signal of a blank sample. Disable by setting the parameter to NULL.

-   *min_feature:* Default set to 6, according to EMA guidelines 2012. Minimal number of consecutive signals present per dilution/concentration curve, to consider this curve as linear.

-   *output_name*: name of the output file

-   *output_dir: directory where the output files should be saved*

### With these arguments you can run the example data set:

``` r
library(MSlineaR)

# load example data
data_tbl <- MSlineaR::targeted_MS


targetedMSCal <- AssessLinearity(
  analysis_type ="targeted",
  input_data = data_tbl,
  column_sample_type = "Sample.Type" ,
  sample_type_QC = c("Pooled QC","Reference QC"),
  sample_type_sample = "Sample",
  sample_type_serial = "Calibration Standard",
  sample_type_blank = "Blank",
  column_sample_ID = "Sample.Identification",

  column_compound_ID = "Compound",
  column_batch = "Batch",
  column_X = "Concentration",
  dilution_factor = 3,
  signal_blank_ratio = 5,
  column_Y = "Area",
  column_Y_sample = "Area",
  min_feature = 5,
  output_name = "Test_targeted",
  output_dir = "S:/Projects/2022/MSTARS/Janine_linearity_project/Benchmark data/targeted/230419",
  nCORE = 4

)
```

## Results

## Optional Input Arguments:

-   *sample_type_QC :*

-   *column_Batch = "Batch"*

-   *column_Y\_sample = "Area"*

-   *column_class_sample = NULL*

-   *#linear_range*

-   *nCORE = 4*

-   *get_output = TRUE*

-   *#transform_data*

-   *transform = c(TRUE, FALSE)[1]*

    -   *transform_x = "log"*

    -   *inverse_x = "exp"*

    -   *transform_y = "log"*

    -   *inverse_y = "exp"*

-   *#first_outlier_detection*

-   *first_outlier_detection = c(TRUE, FALSE)[1]*

    -   *FOD_model = c("logistic", "linear", "quadratic")*

    -   *FOD_sdres_min = 1*

    -   *FOD_stdres_max = 2*

-   *trimming = c(TRUE, FALSE)[1]*

-   *#second_outlier_detection*

-   *second_outlier_detection = c(TRUE, FALSE)[1]*

    -   *SOD_model = c("logistic", "linear", "quadratic")*

    -   *SOD_sdres_min = 1*

    -   *SOD_stdres_max = 2*

    -   *R2_min = 0.9*

-   *#linear_range*

-   *LR_sd_res_factor = 2*

-   *Batch_harmonization = c(TRUE, FALSE)[1]*

-   *calculate_concentration = c(TRUE, FALSE)[1]*

-   *get_linearity_status_samples = c(TRUE, FALSE)[1]*

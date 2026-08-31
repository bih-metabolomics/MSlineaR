test_that("MS_AssessLinearity works with dilution series, Samples and QC only ", {

  sample_test <- MSlineaR::Sample_tbl |>
    dplyr::filter(
      Sample.Type %in% c("Calibration Standard", "Pooled QC", "Sample")
    )

  # use only one batch
  test_batch <- unique(sample_test$Batch)[1]

  sample_test <- sample_test |>
    dplyr::filter(
      Batch == test_batch
    )

  # use only one feature
  test_feature <- unique(
    MSlineaR::Feature_tbl_long$Compound
  )[1:5]

  feature_test <- MSlineaR::Feature_tbl_long |>
    dplyr::filter(
      Compound %in% test_feature,
      Sample.Identification %in%
        sample_test$Sample.Identification,
      Batch == test_batch
    )

  result <- MS_AssessLinearity(
    analysisType = "untargeted",

    inputData_feature = feature_test,
    inputData_sample = sample_test,

    column_sampleType = "Sample.Type",

    sampleType_serial = "Calibration Standard",
    sampleType_sample = "Sample",
    sampleType_QC = "Pooled QC",
    sampleType_blank = NULL,
    signal_blank_ratio = NULL,

    column_sampleID = "Sample.Identification",
    column_featureID = "Compound",
    column_injectionOrder = "Sequence.Position",
    column_batch = "Batch",
    column_X = "Dilution",
    column_Y = "Area",
    column_sampleClass = "Group",

    min_feature = 5,
    nCORE = 1,

    output_name = "test_dilution_only",
    output_dir = tempdir(),

    get_output = TRUE,
    which_output = "DilutionCurves"
  )

  expect_type(result, "list")

  expect_named(
    result,
    c(
      "All_DilutionCurves_Signals",
      "All_DilutionCurves_Features",
      "result"
    )
  )

  expect_s3_class(
    result$All_DilutionCurves_Signals,
    "data.frame"
  )

  expect_s3_class(
    result$All_DilutionCurves_Features,
    "data.frame"
  )

  expect_s3_class(
    result$result,
    "data.frame"
  )

  expect_gt(
    nrow(result$All_DilutionCurves_Signals),
    0
  )

  expect_gt(
    nrow(result$All_DilutionCurves_Features),
    0
  )

  expect_gt(
    nrow(result$result),
    0
  )
})

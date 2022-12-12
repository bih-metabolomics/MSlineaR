test_that("checks chooseModel",{
  dat <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                           ID = rep(1,10),
                                           REPLICATE = rep(1,10),
                                           X = 1:10,
                                           Y = 1:10,
                                           DilutionPoint = 1:10))

  dat2 <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                           ID = rep(1,10),
                                           REPLICATE = rep(1,10),
                                           X = 1:10,
                                           Y = c(1:5,4:0),
                                           DilutionPoint = 1:10))



  y = "Y"
  x = "X"

  testthat::expect_equal(chooseModel(dat, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "linear")
  testthat::expect_equal(chooseModel(dat, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$cor, 1)

  testthat::expect_equal(chooseModel(dat2, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$model.name, "quadratic")
  testthat::expect_equal(round(chooseModel(dat2, tidyselect::all_of(y), tidyselect::all_of(x), model)[[1]]$cor, 3), 0.962)


})



test_that("checks outlier", {

  dat <- data.table::data.table(data.frame(groupIndices = rep(1,10),
                                           ID = rep(1,10),
                                           REPLICATE = rep(1,10)
                                           X = 1:10,
                                           Y = 1:10,
                                           DilutionPoint = 1:10))


})

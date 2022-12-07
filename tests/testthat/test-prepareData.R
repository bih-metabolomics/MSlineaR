testthat::test_that("function prepareData works", {

  dat = data.table::data.table(data.frame(ID = 1:10, X = 1:10, Y = 1:10))
  dat2 = data.table::data.table(data.frame(ID = 1:10, X = -1:8, Y = 1:10))
  dat3 = data.table::data.table(data.frame(ID = 1:10, X = 1:10, Y = -1:8))

  COLNAMES = c(ID = "ID", REPLICATE = NULL, X = "X", Y = "Y")
  nCORE = 1
  LOG_TRANSFORM = TRUE
  MIN_FEATURE = 3
  testthat::expect_vector(prepareData(dat)$ID, ptype = character(), size= 10)
  testthat::expect_vector(prepareData(dat)$REPLICATE, ptype = character(), size= 10)
  testthat::expect_vector(prepareData(dat)$X, ptype = double(), size= 10)
  testthat::expect_vector(prepareData(dat)$Y, ptype = double(), size= 10)
  testthat::expect_error(prepareData(dat2))
  testthat::expect_error(prepareData(dat3))
  testthat::expect_error(prepareData(dat, MIN_FEATURE = 2))
  testthat::expect_error(prepareData(dat, nCORE = -1))
  testthat::expect_error(prepareData(dat, LOG_TRANSFORM = "TRUE"))

  })





